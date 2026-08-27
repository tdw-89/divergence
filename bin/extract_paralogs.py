#!/usr/bin/env python
"""Select orthogroups containing paralogs of one or more target species.

Unlike the original version of this script, the emitted FASTA contains the
*whole* orthogroup (every species), not just the target species' genes. The
extra sequences are not reported on: they exist to break up long branches and
to stabilise branch-length estimation in the downstream codon model. The genes
we actually want pairwise dS for are listed separately in a manifest file.

For each selected orthogroup three files are written:
  <OG>.faa           all orthogroup sequences (unaligned protein)
  <OG>.tree          the OrthoFinder gene tree, if one exists
  <OG>.paralogs.txt  target-species gene IDs, as `species<TAB>gene`

plus a single `orthogroup_summary.tsv` describing every orthogroup considered.
"""

import argparse
import csv
import os
import re
import sys
from functools import lru_cache
from io import StringIO

from Bio import Phylo, SeqIO

SUMMARY_FILE = "orthogroup_summary.tsv"

# OrthoFinder 3.x writes every resolved gene tree into one file, one per line,
# as "OG0000000: (newick);". Older versions wrote a file per orthogroup, so both
# layouts are supported.
COMBINED_TREE_FILES = ("Resolved_Gene_Trees.txt", "Gene_Trees.txt")
TREE_PATTERNS = ("{og}_tree.txt", "{og}_tree.txt.rooted", "{og}.txt", "{og}.tree")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select orthogroups containing paralogs of one or more target species."
    )
    parser.add_argument("--tsv", required=True, help="Path to Orthogroups.tsv")
    parser.add_argument(
        "--seqdir", required=True, help="Path to the Orthogroup_Sequences directory"
    )
    parser.add_argument(
        "--treedir",
        help="Path to the Resolved_Gene_Trees directory. Orthogroups with no tree "
        "are still emitted; the downstream step falls back to pairwise-only estimation.",
    )
    parser.add_argument(
        "--species",
        required=True,
        help="Comma-separated target species names, each matched as a prefix against "
        "the TSV header, or 'all' for every species in the table",
    )
    parser.add_argument(
        "--min-paralogs",
        type=int,
        default=2,
        help="Minimum genes a target species needs in an orthogroup to be reported "
        "on (default: 2). An orthogroup is kept if any target species reaches it.",
    )
    parser.add_argument(
        "--max-seqs",
        type=int,
        default=0,
        help="Skip orthogroups with more than this many sequences in total. "
        "0 (the default) means no limit. Guards against gene families whose "
        "pairwise cost would dominate a batch.",
    )
    parser.add_argument("--outdir", default=".", help="Directory to write output files into")
    return parser.parse_args()


def find_species_column(header, species):
    """Return the index of the header column naming `species`.

    An exact match wins outright. Species are named by the samplesheet's
    `species` column and the proteome is written under that name, so the header
    spells it exactly; preferring the exact hit stops a name that happens to
    prefix another -- `human` beside `human_alt` -- from reading as ambiguous
    and failing a run whose intent was never in doubt.

    Prefix matching stays as the fallback because `--orthofinder_dir` can supply
    a run whose species were named from filenames, where the header carries
    whatever suffix the file had (`GCF_049306965.2_protein` for
    `GCF_049306965.2`). It still errors on zero or multiple matches.
    """
    exact = [i for i, col in enumerate(header) if col == species]
    if len(exact) == 1:
        return exact[0]
    matches = exact or [i for i, col in enumerate(header) if col.startswith(species)]
    if not matches:
        raise ValueError(
            f"Species '{species}' does not name a column in Orthogroups.tsv. "
            f"Species are named by the samplesheet's `species` column, and this "
            f"run has: {', '.join(header[1:])}"
        )
    if len(matches) > 1:
        raise ValueError(
            f"Species '{species}' matches multiple columns: {[header[i] for i in matches]}"
        )
    return matches[0]


def find_species_columns(header, species_arg):
    """Resolve a comma-separated species list (or 'all') to {column name: index}.

    Each name is matched as a prefix and must match exactly one column, so a
    typo or an ambiguous prefix is an error rather than a silently empty result.
    """
    if species_arg.strip().lower() == "all":
        return {header[i]: i for i in range(1, len(header))}

    requested = [name.strip() for name in species_arg.split(",") if name.strip()]
    if not requested:
        raise ValueError("No target species given")

    resolved = {}
    for name in requested:
        index = find_species_column(header, name)
        if header[index] in resolved:
            raise ValueError(f"Species '{name}' resolves to an already-selected column")
        resolved[header[index]] = index
    return resolved


def parse_genes(cell):
    """OrthoFinder separates genes within a cell with a comma and a space."""
    cell = cell.strip()
    if not cell:
        return []
    return [gene for gene in (g.strip() for g in cell.split(",")) if gene]


def load_combined_trees(treedir):
    """Read OrthoFinder's one-file-per-run tree output into {orthogroup: newick}."""
    trees = {}
    if not treedir:
        return trees
    for name in COMBINED_TREE_FILES:
        path = os.path.join(treedir, name)
        if not os.path.exists(path):
            continue
        with open(path) as handle:
            for line in handle:
                line = line.strip()
                if not line:
                    continue
                # Orthogroup ids contain no colon, so the first one splits the line.
                og_id, sep, newick = line.partition(":")
                if sep and newick.strip():
                    trees[og_id.strip()] = newick.strip()
        break
    return trees


def locate_tree(treedir, og_id):
    """Fall back to the per-orthogroup tree files older OrthoFinder versions wrote."""
    if not treedir:
        return None
    for pattern in TREE_PATTERNS:
        candidate = os.path.join(treedir, pattern.format(og=og_id))
        if os.path.exists(candidate):
            return candidate
    return None


@lru_cache(maxsize=None)
def prefix_pattern(species):
    """A regex matching `<species>_` in any spelling OrthoFinder may have used.

    OrthoFinder writes the species name verbatim in the Orthogroups.tsv header
    but rewrites some punctuation to underscores in gene-tree tip labels, so a
    species whose name contains dots -- an accession such as
    `GCF_049306965.2_GRCz12tu` -- is spelt differently in the two files.

    Which characters get rewritten is not worth guessing: observed output
    replaces `.` but leaves `-` alone (`GCF_053564925.1_Olatipes_Hd-rR_3.1`
    becomes `GCF_053564925_1_Olatipes_Hd-rR_3_1`). So instead of hard-coding a
    character class, every non-alphanumeric character is allowed to match
    either itself or an underscore. Over-matching is harmless: the caller only
    strips a prefix when the remainder is a gene genuinely in the orthogroup.

    The gene id that follows is *not* rewritten, so only the species side is
    treated this way.
    """
    parts = [
        re.escape(ch) if ch.isalnum() or ch == "_" else f"[{re.escape(ch)}_]"
        for ch in species
    ]
    return re.compile("".join(parts) + "_")


def normalise_tree(newick, genes, species_names):
    """Strip OrthoFinder's species prefixes from tip labels.

    Gene trees label tips `<species>_<gene id>` while Orthogroups.tsv and the
    orthogroup FASTAs use the bare gene id, so the labels have to be brought
    back into line before the tree can be matched to an alignment. A prefix is
    only removed when what remains is a gene actually in this orthogroup, which
    keeps the rewrite from mangling ids that happen to look like a prefix.

    Returns (tree, unmatched_tip_labels).
    """
    tree = Phylo.read(StringIO(newick), "newick")
    patterns = [prefix_pattern(species) for species in species_names]
    unmatched = []
    for tip in tree.get_terminals():
        name = tip.name
        if name is None:
            continue
        if name in genes:
            continue
        for pattern in patterns:
            match = pattern.match(name)
            if match and name[match.end():] in genes:
                tip.name = name[match.end():]
                break
        else:
            unmatched.append(name)
    return tree, unmatched


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    summary_rows = []
    written = 0

    combined_trees = load_combined_trees(args.treedir)

    with open(args.tsv, "r") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        target_columns = find_species_columns(header, args.species)
        # Every column after the orthogroup ID is a species.
        species_columns = list(range(1, len(header)))
        species_names = [header[i] for i in species_columns]

        for row in reader:
            og_id = row[0]

            # Genes per species, read once: a row can be short if OrthoFinder
            # omitted trailing empty columns, so every cell is read through the
            # same bounds check rather than only the target ones.
            genes_by_column = {
                i: parse_genes(row[i]) if i < len(row) else [] for i in species_columns
            }

            target_genes = {
                name: genes_by_column[index]
                for name, index in target_columns.items()
                if len(genes_by_column[index]) >= args.min_paralogs
            }

            n_total = sum(len(g) for g in genes_by_column.values())
            n_species = sum(1 for g in genes_by_column.values() if g)
            n_target_paralogs = sum(len(g) for g in target_genes.values())

            status = "kept"
            if not target_genes:
                status = "too_few_paralogs"
            elif args.max_seqs and n_total > args.max_seqs:
                status = "too_many_sequences"

            if status == "kept":
                seq_file = os.path.join(args.seqdir, f"{og_id}.fa")
                if not os.path.exists(seq_file):
                    print(
                        f"WARNING: Sequence file {seq_file} not found. Skipping {og_id}.",
                        file=sys.stderr,
                    )
                    status = "missing_sequences"

            wrote_tree = False
            # Why the tree is (or is not) there. `has_tree = no` on its own
            # cannot tell an orthogroup OrthoFinder never built a tree for --
            # it writes none below four sequences, which is most small
            # orthogroups and entirely benign -- from one whose tree we had and
            # then threw away. The latter has silently degraded whole runs to
            # pairwise-only twice, so it is recorded rather than inferred.
            tree_status = "not_applicable"

            if status == "kept":
                # Emit the whole orthogroup, every species included.
                records = list(SeqIO.parse(seq_file, "fasta"))
                if len(records) < 2:
                    status = "too_few_sequences"
                else:
                    out_fasta = os.path.join(args.outdir, f"{og_id}.faa")
                    SeqIO.write(records, out_fasta, "fasta")

                    manifest = os.path.join(args.outdir, f"{og_id}.paralogs.txt")
                    with open(manifest, "w") as fh:
                        fh.write("species\tgene\n")
                        for name, genes in target_genes.items():
                            for gene in genes:
                                fh.write(f"{name}\t{gene}\n")

                    # Always emit a tree file, empty when no usable tree exists.
                    # Keeping the file set uniform lets the workflow join
                    # alignments, trees and manifests without special-casing
                    # absent trees, and an empty file tells the dS step to fall
                    # back to pairwise-only estimation.
                    out_tree = os.path.join(args.outdir, f"{og_id}.tree")
                    newick = combined_trees.get(og_id)
                    if newick is None:
                        per_og = locate_tree(args.treedir, og_id)
                        if per_og:
                            newick = open(per_og).read().strip()

                    if not newick:
                        tree_status = "none_from_orthofinder"
                    else:
                        genes_in_og = {
                            gene for genes in genes_by_column.values() for gene in genes
                        }
                        try:
                            tree, unmatched = normalise_tree(newick, genes_in_og, species_names)
                        except Exception as exc:
                            print(
                                f"WARNING: Could not parse the gene tree for {og_id} ({exc}).",
                                file=sys.stderr,
                            )
                            tree, unmatched = None, []
                            tree_status = "unparseable"
                        if tree is not None and unmatched:
                            print(
                                f"WARNING: {len(unmatched)} tip(s) in the {og_id} gene tree do "
                                f"not correspond to any gene in the orthogroup "
                                f"(e.g. {unmatched[:2]}); dropping the tree.",
                                file=sys.stderr,
                            )
                            tree = None
                            tree_status = "unmatched_tips"
                        if tree is not None:
                            Phylo.write(tree, out_tree, "newick")
                            wrote_tree = True
                            tree_status = "ok"

                    if not wrote_tree:
                        open(out_tree, "w").close()

                    written += 1

            summary_rows.append(
                {
                    "orthogroup": og_id,
                    "n_sequences": n_total,
                    "n_species": n_species,
                    "n_target_species": len(target_genes),
                    "n_target_paralogs": n_target_paralogs,
                    "target_paralog_counts": ";".join(
                        f"{name}:{len(genes)}" for name, genes in sorted(target_genes.items())
                    ),
                    "has_tree": "yes" if wrote_tree else "no",
                    "tree_status": tree_status,
                    "status": status,
                }
            )

    summary_path = os.path.join(args.outdir, SUMMARY_FILE)
    with open(summary_path, "w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            lineterminator="\n",
            fieldnames=[
                "orthogroup",
                "n_sequences",
                "n_species",
                "n_target_species",
                "n_target_paralogs",
                "target_paralog_counts",
                "has_tree",
                "tree_status",
                "status",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    print(
        f"Selected {written} orthogroup(s) in which at least one target species has "
        f">={args.min_paralogs} genes, out of {len(summary_rows)} considered.",
        file=sys.stderr,
    )

    if written == 0:
        print(
            f"ERROR: No orthogroup contained at least {args.min_paralogs} genes for any "
            f"of the target species ({args.species}). Check that --target_species names "
            f"columns in Orthogroups.tsv.",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
