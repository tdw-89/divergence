#!/usr/bin/env python3

"""Estimate pairwise dS (Ks) between the paralogs of one or more target species.

The input alignment holds an entire orthogroup -- every species, not just the
targets. Pairs are formed within each target species, so only paralogs are
reported, never orthologs. Sequences from non-target species still take part in
the alignment and the model fit: they break up long branches so that the codon
model's multiple-hit correction has something to work with.

The fitted tree is written out as `<OG>_dS.nwk` and `<OG>_dN.nwk`, with branch
lengths in dS and dN units and the original sequence IDs restored. Any distance
the TSV does not report -- between orthologs, or to an internal node -- can be
read off those.

Three estimates are produced for each paralog pair:

  tree_*  CODEML M0 (`runmode = 0`) fitted to the whole orthogroup on a fixed
          topology. dS is the sum of the per-branch dS values CODEML reports
          along the path connecting the two paralogs. This is the primary
          estimate: each branch is individually short, and its length is
          informed by every sequence in the family.
  pair_*  CODEML pairwise ML (`runmode = -2`) on just the two sequences, after
          dropping codon columns gapped in either one. Independent of the tree,
          so disagreement with tree_* flags saturation, a bad alignment or a
          bad topology.
  yn00_*  Yang & Nielsen (2000) counting method on the same two-sequence
          alignment. A different method class again, so it catches problems
          both ML estimates share.

No pair is ever dropped for being low quality. Quality columns are emitted
alongside the estimates and curation is left to downstream analysis.

CODEML is driven through a control file and its output is parsed here rather
than through Biopython's PAML wrapper, whose parser fails on truncated result
tables. Biopython is still used for sequence, alignment and tree handling.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import random
import re
import shutil
import subprocess
import sys
import tempfile
import warnings
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from io import StringIO
from pathlib import Path

from Bio import AlignIO, Phylo, SeqIO
from Bio.Align import MultipleSeqAlignment
from Bio.Data import CodonTable
from Bio.SeqRecord import SeqRecord

# Bio.codonalign warns on import that it is experimental. We know.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio.codonalign import build as build_codon_alignment

# PAML's `icode` numbering is its own; map it to the NCBI translation table
# ids Biopython uses so the codon alignment is built under the same code
# CODEML will assume. PAML icode 3 is the Mycoplasma/Spiroplasma code.
PAML_ICODE_TO_NCBI = {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 9, 7: 10, 8: 12, 9: 13, 10: 15}

GAP = "-"

# OrthoFinder writes no gene tree below this many sequences, and `load_and_prepare_tree`
# would reject one anyway. Below it there is no tree-based estimate to be had for
# reasons that have nothing to do with the data being difficult.
MIN_TREE_SEQS = 4

_NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?|nan|-?inf"
_RE_PAIRWISE = re.compile(
    r"t=\s*(?P<t>" + _NUM + r")\s+"
    r"S=\s*(?P<S>" + _NUM + r")\s+"
    r"N=\s*(?P<N>" + _NUM + r")\s+"
    r"dN/dS=\s*(?P<omega>" + _NUM + r")\s+"
    r"dN\s*=\s*(?P<dN>" + _NUM + r")\s+"
    r"dS\s*=\s*(?P<dS>" + _NUM + r")"
)
# CODEML's per-branch table, printed under "dN & dS for each branch":
#    branch          t       N       S   dN/dS      dN      dS  N*dN  S*dS
#      7..8      0.021   880.2   319.8  0.0462  0.0011  0.0228   0.9   7.3
_RE_BRANCH = re.compile(
    r"^\s*(?P<parent>\d+)\.\.(?P<child>\d+)\s+"
    r"(?P<t>" + _NUM + r")\s+"
    r"(?P<N>" + _NUM + r")\s+"
    r"(?P<S>" + _NUM + r")\s+"
    r"(?P<omega>" + _NUM + r")\s+"
    r"(?P<dN>" + _NUM + r")\s+"
    r"(?P<dS>" + _NUM + r")\s+"
    r"(?P<NdN>" + _NUM + r")\s+"
    r"(?P<SdS>" + _NUM + r")\s*$",
    re.M,
)

# YN00's result row, under the "(B) Yang & Nielsen (2000) method" heading:
# seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE
#    2    1   370.8   829.2   0.0694  4.6000  0.2265 0.0112 +- 0.0037  0.0497 +- 0.0119
_RE_YN00 = re.compile(
    r"^\s*\d+\s+\d+\s+"
    r"(?P<S>" + _NUM + r")\s+"
    r"(?P<N>" + _NUM + r")\s+"
    r"(?P<t>" + _NUM + r")\s+"
    r"(?P<kappa>" + _NUM + r")\s+"
    r"(?P<omega>" + _NUM + r")\s+"
    r"(?P<dN>" + _NUM + r")\s*\+-\s*(?P<dN_SE>" + _NUM + r")\s+"
    r"(?P<dS>" + _NUM + r")\s*\+-\s*(?P<dS_SE>" + _NUM + r")",
    re.M,
)
_YN00_SECTION = "(B) Yang & Nielsen (2000) method"
_BRANCH_SECTION = "dN & dS for each branch"

_RE_OMEGA = re.compile(r"omega\s*\(dN/dS\)\s*=\s*(" + _NUM + r")")
_RE_KAPPA = re.compile(r"kappa\s*\(ts/tv\)\s*=\s*(" + _NUM + r")")
_RE_LNL = re.compile(r"lnL\s*\([^)]*\)\s*:\s*(" + _NUM + r")")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Estimate pairwise dS between the paralogs of one or more target species with CODEML."
    )
    parser.add_argument("-m", "--msa", required=True, help="Aligned protein FASTA for the orthogroup.")
    parser.add_argument(
        "-c",
        "--cds",
        required=True,
        nargs="+",
        help="Nucleotide CDS FASTA(s). One per species; IDs must match the protein IDs exactly.",
    )
    parser.add_argument(
        "-p",
        "--paralogs",
        required=True,
        help="Manifest of target genes as `species<TAB>gene` rows. Pairs are formed "
        "within each species, so only paralogs are reported.",
    )
    parser.add_argument("-o", "--output", required=True, help="Output TSV path.")
    parser.add_argument(
        "-t",
        "--tree",
        help="Newick gene tree for the orthogroup. Without it the tree-based "
        "estimate is skipped and only the pairwise columns are filled in.",
    )
    parser.add_argument(
        "--orthogroup",
        help="Orthogroup identifier for the output. Defaults to the MSA filename stem.",
    )
    parser.add_argument(
        "--icode",
        type=int,
        default=0,
        help="PAML genetic code index (default: 0, the universal code). "
        "3 is the Mycoplasma/Spiroplasma code.",
    )
    parser.add_argument("--msa-format", default="fasta", help="Format of --msa (default: fasta).")
    parser.add_argument("--codeml-command", default="codeml", help="CODEML executable (default: codeml).")
    parser.add_argument("--yn00-command", default="yn00", help="YN00 executable (default: yn00).")
    parser.add_argument(
        "--skip-tree", action="store_true", help="Skip the tree-based estimate even if a tree is given."
    )
    parser.add_argument("--skip-pairwise", action="store_true", help="Skip the pairwise CODEML estimate.")
    parser.add_argument("--skip-yn00", action="store_true", help="Skip the YN00 estimate.")
    parser.add_argument(
        "--codeml-method",
        type=int,
        default=0,
        help="CODEML `method`: 0 simultaneous, 1 one branch at a time (faster on big trees).",
    )
    parser.add_argument(
        "--codeml-timeout",
        type=int,
        default=0,
        help="Abandon a single CODEML or YN00 invocation after this many seconds "
        "(0, the default, waits indefinitely). Deeply saturated sequences make the "
        "ML optimiser crawl toward an unbounded dS, so one orthogroup can otherwise "
        "consume an entire task's time budget. A timed-out estimate is reported as "
        "NA and the rest of the orthogroup still runs.",
    )
    parser.add_argument(
        "--max-pairwise-pairs",
        type=int,
        default=0,
        help="Run the pairwise CODEML and YN00 cross-checks on at most this many "
        "pairs per orthogroup, chosen at random (0, the default, runs all). The "
        "tree-based dS -- the primary estimate -- comes from a single M0 fit and "
        "covers every pair regardless, so this caps a cost that grows with the "
        "square of the paralog count while the estimate it guards does not. Lifted "
        "only for orthogroups too small for a gene tree to exist, where the pairwise "
        "columns are the only result there is and there are few pairs anyway; a tree "
        "that was available and failed to fit stays capped. Sampling is seeded on the "
        "orthogroup name, so a rerun picks the same pairs.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Run this many pairwise cross-checks concurrently (default: 1). The "
        "work is external CODEML and YN00 processes, so threads are enough.",
    )
    parser.add_argument(
        "--keep-temp", action="store_true", help="Keep CODEML working directories for debugging."
    )
    return parser.parse_args()


# --------------------------------------------------------------------------
# Input handling
# --------------------------------------------------------------------------


def read_protein_alignment(msa_path: Path, msa_format: str) -> MultipleSeqAlignment:
    alignment = AlignIO.read(str(msa_path), msa_format)
    seen: set[str] = set()
    duplicates: set[str] = set()
    for record in alignment:
        if record.id in seen:
            duplicates.add(record.id)
        seen.add(record.id)
    if duplicates:
        raise ValueError(
            "Duplicate sequence ID(s) in the protein alignment: " + ", ".join(sorted(duplicates))
        )
    return alignment


class CdsIndex(Mapping):
    """Lazy, offset-based lookup of CDS records across several files.

    An orthogroup needs a few dozen of the hundreds of thousands of records in
    the run's CDS files, so parsing them all into memory costs half a gigabyte
    per process to answer a handful of questions -- and CODEML_BATCH runs
    several `ks.py` processes concurrently, each paying it again. `SeqIO.index`
    stores only a file offset per record and parses on lookup.

    Presents the same read-only mapping interface the eager dict did, so
    `.get()`, `in` and iteration over ids all still work.
    """

    __slots__ = ("_indexes",)

    def __init__(self, indexes: list[Mapping]) -> None:
        self._indexes = indexes

    def __getitem__(self, key: str) -> SeqRecord:
        for index in self._indexes:
            if key in index:
                return index[key]
        raise KeyError(key)

    def __iter__(self):
        for index in self._indexes:
            yield from index

    def __len__(self) -> int:
        return sum(len(index) for index in self._indexes)

    def close(self) -> None:
        """Release the open file handles the offset indexes hold."""
        for index in self._indexes:
            index.close()


def build_cds_index(cds_paths: list[Path]) -> CdsIndex:
    """Index every CDS file by record id, exactly as written.

    Matching is deliberately strict: an earlier attempt at fuzzy ID parsing was
    reverted because it silently paired the wrong transcripts. A duplicate ID
    across two species' files is an error rather than a silent overwrite.
    """
    indexes: list[Mapping] = []
    origin: dict[str, Path] = {}
    try:
        for path in cds_paths:
            # Raises on a duplicate id *within* one file; across files is our job.
            index = SeqIO.index(str(path), "fasta")
            indexes.append(index)
            for key in index:
                if key in origin:
                    raise ValueError(
                        f"Duplicate CDS ID '{key}' found in both {origin[key]} "
                        f"and {path}. IDs must be unique across all CDS files."
                    )
                origin[key] = path
        if not origin:
            raise ValueError(
                f"No sequences read from CDS file(s): {', '.join(str(p) for p in cds_paths)}"
            )
    except Exception:
        # Don't leave the handles of the files we did open dangling.
        for index in indexes:
            index.close()
        raise
    return CdsIndex(indexes)


def read_manifest(path: Path) -> dict[str, list[str]]:
    """Read `species<TAB>gene` rows into {species: [gene, ...]}.

    Grouping by species is what keeps cross-species (ortholog) pairs out of the
    output: pairs are only ever formed within one species' gene list.
    """
    genes_by_species: dict[str, list[str]] = {}
    first = True
    for lineno, raw in enumerate(path.read_text().splitlines(), start=1):
        line = raw.strip()
        if not line:
            continue
        parts = [part.strip() for part in line.split("\t")]
        if len(parts) != 2 or not all(parts):
            raise ValueError(
                f"{path}:{lineno}: expected 'species<TAB>gene', got {raw!r}"
            )
        if first and parts == ["species", "gene"]:
            first = False
            continue
        first = False
        genes_by_species.setdefault(parts[0], []).append(parts[1])
    if not genes_by_species:
        raise ValueError(f"Paralog manifest {path} is empty.")
    return genes_by_species


def match_cds(protein_ids: list[str], index: Mapping) -> tuple[list[SeqRecord], list[str]]:
    matched: list[SeqRecord] = []
    missing: list[str] = []
    for protein_id in protein_ids:
        record = index.get(protein_id)
        if record is None:
            missing.append(protein_id)
        else:
            matched.append(record)
    return matched, missing


# --------------------------------------------------------------------------
# Alignment helpers
# --------------------------------------------------------------------------


def select_reportable(
    paralogs: dict[str, list[str]], retained: set[str]
) -> dict[str, list[str]]:
    """Target-species gene lists restricted to sequences we can still analyse.

    Pairs are formed within each species' list, which is what keeps ortholog
    pairs out of the output, so a species with fewer than two surviving genes
    contributes nothing and is dropped.
    """
    reportable: dict[str, list[str]] = {}
    for species, genes in paralogs.items():
        kept = [gene for gene in genes if gene in retained]
        if len(kept) >= 2:
            reportable[species] = kept
    return reportable


def crosscheck_cap_applies(n_seqs: int, n_pairs: int, cap: int) -> bool:
    """Whether `--max-pairwise-pairs` should be enforced for this orthogroup.

    The M0 fit gives dS for every pair from one optimisation, so the primary
    estimate costs the same whether a family has ten paralogs or four hundred.
    The pairwise and YN00 cross-checks do not: they are two external processes
    per pair, and pair count grows quadratically. A 417-paralog family is 86,751
    pairs, roughly an hour and a half of subprocess spawning, against minutes for
    the fit it is checking. Their job is to say whether the estimate sits in a
    trustworthy regime, and a few thousand pairs establish that as well as all of
    them, so the count is capped.

    The cap is lifted only for orthogroups too small for a tree to have existed
    at all -- OrthoFinder writes none below MIN_TREE_SEQS sequences -- where the
    pairwise columns are the only result there is and there are at most three
    pairs to run. It is deliberately *not* lifted when a tree was available and
    the fit failed (a timeout, MAXNSONS, tips that would not reconcile). Those
    are the large families, so keying on whether the fit *succeeded* would run
    the full quadratic cross-check on exactly the orthogroups that can least
    afford it -- and `manifestPairCost`, which packs the batches, has already
    predicted the capped cost for them.
    """
    if not cap or n_pairs <= cap:
        return False
    return n_seqs >= MIN_TREE_SEQS


def build_codon_alignment_tolerantly(
    alignment: MultipleSeqAlignment, cds_records: list[SeqRecord], codon_table
) -> tuple[MultipleSeqAlignment, list[str]]:
    """Build the codon alignment, dropping only the sequences that will not reconcile.

    `Bio.codonalign.build` is all or nothing: one protein its CDS does not
    translate back to -- an annotated frameshift, a partial or low-quality CDS,
    a selenocysteine -- raises and takes the entire orthogroup with it. RefSeq
    carries such records at a percent or two of sequences, and at 1-(1-p)^n that
    silently destroys most of the *large* orthogroups, which are exactly the
    ones worth having.

    The whole-orthogroup build is tried first, so a clean orthogroup costs
    nothing extra. Only when that fails is each sequence screened on its own,
    using Biopython's own build as the test rather than reimplementing its
    notion of a match. Returns (codon_alignment, dropped_ids).
    """

    def attempt(aln: MultipleSeqAlignment, records: list[SeqRecord]):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return build_codon_alignment(aln, records, codon_table=codon_table)

    try:
        return attempt(alignment, cds_records), []
    except Exception:
        pass

    cds_by_id = {record.id: record for record in cds_records}
    good: list[SeqRecord] = []
    dropped: list[str] = []
    for record in alignment:
        cds = cds_by_id.get(record.id)
        if cds is None:
            dropped.append(record.id)
            continue
        try:
            attempt(MultipleSeqAlignment([record]), [cds])
        except Exception:
            dropped.append(record.id)
        else:
            good.append(record)

    if len(good) < 2:
        raise ValueError(
            f"only {len(good)} sequence(s) reconcile with their CDS "
            f"({len(dropped)} dropped)"
        )
    return attempt(MultipleSeqAlignment(good), [cds_by_id[r.id] for r in good]), dropped


def make_short_id_map(sequence_ids: list[str]) -> dict[str, str]:
    """PAML truncates long sequence names, so work in S00001-style ids."""
    return {seq_id: f"S{idx:05d}" for idx, seq_id in enumerate(sequence_ids, start=1)}


def strip_pair_gaps(seq_a: str, seq_b: str) -> tuple[str, str]:
    """Drop codon columns gapped in either sequence.

    Done per pair rather than across the orthogroup on purpose: a single
    fragmentary paralog would otherwise delete columns from every other pair
    in the family.
    """
    kept_a: list[str] = []
    kept_b: list[str] = []
    for start in range(0, min(len(seq_a), len(seq_b)) - 2, 3):
        codon_a = seq_a[start : start + 3]
        codon_b = seq_b[start : start + 3]
        if GAP in codon_a or GAP in codon_b:
            continue
        if len(codon_a) < 3 or len(codon_b) < 3:
            continue
        kept_a.append(codon_a)
        kept_b.append(codon_b)
    return "".join(kept_a), "".join(kept_b)


def ungapped_codons(seq: str) -> int:
    """Codons this sequence actually occupies in the codon alignment."""
    return sum(1 for i in range(0, len(seq) - 2, 3) if seq[i] != GAP)


def pct_identity(seq_a: str, seq_b: str) -> tuple[float, float]:
    """Bidirectional identity; gaps in the opposite sequence count as mismatches."""
    matches = sum(1 for ca, cb in zip(seq_a, seq_b) if ca == cb and ca != GAP)
    non_gap_a = sum(1 for c in seq_a if c != GAP)
    non_gap_b = sum(1 for c in seq_b if c != GAP)
    return (
        matches / non_gap_a if non_gap_a else 0.0,
        matches / non_gap_b if non_gap_b else 0.0,
    )


def write_phylip(records: list[tuple[str, str]], path: Path) -> None:
    """Strict sequential PAML-style PHYLIP.

    Written by hand because Biopython's PHYLIP writers either truncate names or
    interleave, both of which PAML mis-reads.
    """
    if not records:
        raise ValueError("Refusing to write an empty alignment")
    seq_len = len(records[0][1])
    for name, seq in records:
        if len(seq) != seq_len:
            raise ValueError("Alignment rows differ in length")
    with path.open("w") as handle:
        handle.write(f" {len(records)} {seq_len}\n")
        for name, seq in records:
            handle.write(f"{name}\n{seq.upper()}\n")


# --------------------------------------------------------------------------
# Tree handling
# --------------------------------------------------------------------------


def load_and_prepare_tree(tree_path: Path, keep_ids: set[str], id_map: dict[str, str]):
    """Prune the gene tree to the retained sequences and relabel to short ids.

    Returns None when the tree cannot be reconciled with the alignment, which
    makes the caller fall back to pairwise-only estimation rather than fail.
    """
    tree = Phylo.read(str(tree_path), "newick")
    tip_names = {tip.name for tip in tree.get_terminals() if tip.name}

    missing_from_tree = keep_ids - tip_names
    if missing_from_tree:
        raise ValueError(
            f"{len(missing_from_tree)} alignment sequence(s) absent from the gene tree "
            f"(e.g. {sorted(missing_from_tree)[:3]})"
        )

    for tip in list(tree.get_terminals()):
        if tip.name not in keep_ids:
            tree.prune(tip)

    remaining = [tip for tip in tree.get_terminals() if tip.name]
    if len(remaining) < 3:
        raise ValueError("Fewer than 3 tips remain after pruning; tree adds nothing")

    for tip in remaining:
        tip.name = id_map[tip.name]

    # PAML treats the tree as unrooted under a reversible model; feeding it a
    # rooted (bifurcating) root creates an unidentifiable extra branch.
    _unroot(tree)

    # CODEML cannot read a node with more than three daughters, and OrthoFinder
    # emits them freely, so anything wider has to be split before it is written.
    _resolve_polytomies(tree)

    # OrthoFinder branch lengths are in amino-acid substitution units and make
    # poor starting values for a codon model, so drop them and let CODEML
    # estimate from scratch.
    for clade in tree.find_clades():
        clade.branch_length = None
        if not clade.is_terminal():
            clade.confidence = None
            clade.name = None

    return tree


MAX_PAML_DAUGHTERS = 3


def _resolve_polytomies(tree, root_degree: int = MAX_PAML_DAUGHTERS) -> int:
    """Split multifurcations into binary nodes joined by zero-length branches.

    CODEML has a compile-time cap of three daughter nodes (`MAXNSONS`) and
    aborts with `too many daughter nodes` on anything wider -- so one polytomy
    anywhere in the tree costs the whole orthogroup its tree-based dS, which is
    the primary estimate. OrthoFinder produces polytomies freely: paralogs that
    differ by nothing have zero-length branches between them, and the resolved
    gene tree collapses those rather than inventing an order.

    Re-imposing an arbitrary order is as defensible as the collapse was. The
    branches being restored had length zero, the inserted ones start at zero,
    and CODEML re-estimates every branch from scratch anyway, so the fitted dS
    along any path through the polytomy is unchanged. The root keeps three
    daughters because PAML wants the tree unrooted.

    Children are paired off *balanced* rather than into a ladder: a wide
    polytomy resolved as a ladder is as deep as it is wide, and Biopython walks
    trees recursively, so a few thousand paralogs under one node would overflow
    the stack on the way back out. Pairing halves the count each round, which
    keeps the added depth logarithmic. Returns the number of nodes split.
    """
    split = 0
    for clade in list(tree.find_clades(order="level")):
        limit = root_degree if clade is tree.root else 2
        if len(clade.clades) <= limit:
            continue
        split += 1
        children = clade.clades
        while len(children) > limit:
            paired = [
                type(clade)(clades=[children[i], children[i + 1]], branch_length=0.0)
                for i in range(0, len(children) - 1, 2)
            ]
            if len(children) % 2:
                paired.append(children[-1])
            children = paired
        clade.clades = children
    return split


def _unroot(tree) -> None:
    root = tree.root
    if len(root.clades) != 2:
        return
    # Collapse whichever child is internal into the root to make a trifurcation.
    for i, child in enumerate(root.clades):
        if not child.is_terminal():
            sibling = root.clades[1 - i]
            root.clades = [sibling] + list(child.clades)
            return


def write_tree(tree, path: Path) -> None:
    """Write a topology-only newick for CODEML.

    All branch lengths are stripped: OrthoFinder's are in amino-acid
    substitution units and make poor starting values for a codon model, so
    CODEML (with fix_blength = 0) estimates them from scratch.
    """
    handle = StringIO()
    Phylo.write(tree, handle, "newick")
    newick = handle.getvalue().strip()
    newick = re.sub(r":[0-9.eE+-]+", "", newick)
    with path.open("w") as fh:
        fh.write(newick + "\n")


def parse_branch_table(text: str) -> list[tuple[int, int, float, float]]:
    """Read CODEML's per-branch dN/dS table into (parent, child, dN, dS) edges.

    Under M0, CODEML reports no `dS tree:` newick -- it prints this table with
    branches labelled by node number instead. Node numbers 1..n are the
    sequences in the order they appear in the alignment file, which is the
    order we wrote them in, so tips can be identified without parsing a tree.
    """
    marker = text.find(_BRANCH_SECTION)
    if marker == -1:
        return []
    edges: list[tuple[int, int, float, float]] = []
    for match in _RE_BRANCH.finditer(text[marker:]):
        edges.append(
            (
                int(match.group("parent")),
                int(match.group("child")),
                _to_float(match.group("dN")),
                _to_float(match.group("dS")),
            )
        )
    return edges


_NEWICK_SPECIAL = set("(),:;[]' \t")


def _quote_newick(name: str) -> str:
    if any(ch in _NEWICK_SPECIAL for ch in name):
        return "'" + name.replace("'", "''") + "'"
    return name


def newick_from_branch_table(
    edges: list[tuple[int, int, float, float]],
    name_by_node: dict[int, str],
    weight: str = "dS",
) -> str | None:
    """Rebuild the fitted tree with branch lengths in dS or dN units.

    CODEML prints no dS or dN newick under M0 -- only the per-branch table -- so
    the tree is reassembled from that table. Tips carry their original sequence
    IDs so the result is usable outside this pipeline; distances the TSV does not
    report (ortholog pairs, internal nodes) can be read off it directly.
    """
    pick = (lambda d_n, d_s: d_s) if weight == "dS" else (lambda d_n, d_s: d_n)
    children: dict[int, list[tuple[int, float]]] = {}
    seen_as_child = set()
    for parent, child, d_n, d_s in edges:
        children.setdefault(parent, []).append((child, pick(d_n, d_s)))
        seen_as_child.add(child)

    roots = [node for node in children if node not in seen_as_child]
    if len(roots) != 1:
        return None

    # Iterative post-order: a ladder-shaped gene tree would overflow the
    # recursion limit on a large orthogroup.
    rendered: dict[int, str] = {}
    stack = [(roots[0], False)]
    while stack:
        node, expanded = stack.pop()
        kids = children.get(node)
        if not kids:
            rendered[node] = _quote_newick(name_by_node.get(node, f"node{node}"))
            continue
        if expanded:
            inner = ",".join(f"{rendered[c]}:{w:.6f}" for c, w in kids)
            rendered[node] = f"({inner})"
        else:
            stack.append((node, True))
            for child, _ in kids:
                stack.append((child, False))
    return rendered[roots[0]] + ";"


def tree_path_distances(
    edges: list[tuple[int, int, float, float]],
    node_by_name: dict[str, int],
    from_names: list[str],
) -> dict[tuple[str, str], tuple[float, float]]:
    """Sum dS and dN along the path between each requested pair of tips.

    This is the quantity we actually want: every branch on the path is short
    enough for the multiple-hit correction to work, and each was estimated
    using the whole family rather than just the two sequences.
    """
    adjacency: dict[int, list[tuple[int, float, float]]] = {}
    for parent, child, d_n, d_s in edges:
        adjacency.setdefault(parent, []).append((child, d_n, d_s))
        adjacency.setdefault(child, []).append((parent, d_n, d_s))

    name_by_node = {node: name for name, node in node_by_name.items()}
    distances: dict[tuple[str, str], tuple[float, float]] = {}

    for name in from_names:
        start = node_by_name.get(name)
        if start is None:
            continue
        stack = [(start, 0.0, 0.0)]
        seen = {start}
        while stack:
            node, dn_acc, ds_acc = stack.pop()
            other = name_by_node.get(node)
            if other is not None and node != start:
                distances[(name, other)] = (ds_acc, dn_acc)
            for neighbour, d_n, d_s in adjacency.get(node, []):
                if neighbour not in seen:
                    seen.add(neighbour)
                    stack.append((neighbour, dn_acc + d_n, ds_acc + d_s))
    return distances


# --------------------------------------------------------------------------
# CODEML
# --------------------------------------------------------------------------


CTL_TEMPLATE = """      seqfile = {seqfile}
     treefile = {treefile}
      outfile = {outfile}
        noisy = 0
      verbose = 0
      runmode = {runmode}
      seqtype = 1
    CodonFreq = 2
        clock = 0
       aaDist = 0
        model = 0
      NSsites = 0
        icode = {icode}
        Mgene = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.4
    fix_alpha = 1
        alpha = 0
       Malpha = 0
        ncatG = 1
        getSE = 0
 RateAncestor = 0
   Small_Diff = .5e-6
     cleandata = 0
  fix_blength = 0
       method = {method}
"""


class PamlRun:
    """A finished PAML invocation: its main output file plus its console log.

    PAML reports fatal problems on stdout and *still* leaves a partial main
    output file behind, so a caller that only reads that file sees a truncated
    result with no reason attached. Keeping the console log alongside is what
    turns "output contained no per-branch table" into "too many daughter nodes".
    """

    __slots__ = ("text", "stdout", "returncode")

    def __init__(self, text: str, stdout: str, returncode: int) -> None:
        self.text = text
        self.stdout = stdout
        self.returncode = returncode

    def diagnostic(self, limit: int = 200) -> str:
        """The most informative thing PAML printed, for an error message."""
        for line in self.stdout.splitlines():
            stripped = line.strip()
            if stripped.lower().startswith("error"):
                return stripped[:limit]
        tail = " ".join(self.stdout.split())[-limit:]
        return tail or f"exit status {self.returncode}, no console output"


def _run_paml(command: list[str], workdir: Path, timeout: int) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(
            command,
            cwd=str(workdir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout or None,
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(
            f"{command[0]} exceeded the {timeout}s limit and was abandoned"
        ) from None


def run_codeml(
    workdir: Path,
    phylip: Path,
    treefile: Path,
    runmode: int,
    icode: int,
    method: int,
    command: str,
    timeout: int = 0,
) -> PamlRun:
    """Run CODEML in `workdir` and return its output file plus its console log."""
    outfile = workdir / "mlc"
    ctl = workdir / "codeml.ctl"
    ctl.write_text(
        CTL_TEMPLATE.format(
            seqfile=phylip.name,
            treefile=treefile.name,
            outfile=outfile.name,
            runmode=runmode,
            icode=icode,
            method=method,
        )
    )
    proc = _run_paml([command, ctl.name], workdir, timeout)
    if not outfile.exists():
        raise RuntimeError(
            f"CODEML produced no output file (exit {proc.returncode}). "
            f"Output: {proc.stdout.strip()[-500:]}"
        )
    return PamlRun(outfile.read_text(), proc.stdout, proc.returncode)


def codeml_tree_estimates(
    codon_rows: list[tuple[str, str]],
    tree,
    report_names: list[str],
    workdir: Path,
    icode: int,
    method: int,
    command: str,
    timeout: int = 0,
) -> tuple[
    dict[tuple[str, str], tuple[float, float]], dict[str, float],
    list[tuple[int, int, float, float]], dict[int, str],
]:
    """Fit M0 on the fixed topology and read pairwise dS/dN off the branch table.

    Returns (distances, model statistics, branch edges, {node number: short id}).
    """
    workdir.mkdir(parents=True, exist_ok=True)
    phylip = workdir / "codon_alignment.phy"
    treefile = workdir / "tree.nwk"
    write_phylip(codon_rows, phylip)
    write_tree(tree, treefile)

    run = run_codeml(workdir, phylip, treefile, 0, icode, method, command, timeout)
    text = run.text

    edges = parse_branch_table(text)
    if not edges:
        raise RuntimeError(
            f"CODEML output contained no per-branch dN/dS table -- {run.diagnostic()}"
        )

    # CODEML numbers sequences 1..n in alignment order, which is the order we
    # wrote them, so tip node numbers follow directly from that order.
    node_by_name = {name: index for index, (name, _) in enumerate(codon_rows, start=1)}
    distances = tree_path_distances(edges, node_by_name, report_names)

    stats: dict[str, float] = {}
    for key, pattern in (("omega", _RE_OMEGA), ("kappa", _RE_KAPPA), ("lnL", _RE_LNL)):
        match = pattern.search(text)
        if match:
            stats[key] = _to_float(match.group(1))
    return distances, stats, edges, {index: name for name, index in node_by_name.items()}


def codeml_pairwise_estimate(
    name_a: str,
    seq_a: str,
    name_b: str,
    seq_b: str,
    workdir: Path,
    icode: int,
    command: str,
    timeout: int = 0,
) -> dict[str, float]:
    workdir.mkdir(parents=True, exist_ok=True)
    phylip = workdir / "pair.phy"
    treefile = workdir / "pair.nwk"
    write_phylip([(name_a, seq_a), (name_b, seq_b)], phylip)
    # runmode = -2 ignores the tree, but CODEML still wants the file to exist.
    treefile.write_text(f"({name_a},{name_b});\n")

    run = run_codeml(workdir, phylip, treefile, -2, icode, 0, command, timeout)
    match = _RE_PAIRWISE.search(run.text)
    if not match:
        raise RuntimeError(
            f"CODEML pairwise output contained no result line -- {run.diagnostic()}"
        )
    return {key: _to_float(value) for key, value in match.groupdict().items()}


def yn00_pairwise_estimate(
    name_a: str,
    seq_a: str,
    name_b: str,
    seq_b: str,
    workdir: Path,
    icode: int,
    command: str,
    timeout: int = 0,
) -> dict[str, float]:
    """YN00 on the same two-sequence alignment, parsed from its own output."""
    workdir.mkdir(parents=True, exist_ok=True)
    phylip = workdir / "pair.phy"
    write_phylip([(name_a, seq_a), (name_b, seq_b)], phylip)
    ctl = workdir / "yn00.ctl"
    ctl.write_text(
        f"      seqfile = {phylip.name}\n"
        f"      outfile = yn\n"
        f"      verbose = 0\n"
        f"        icode = {icode}\n"
        f"    weighting = 0\n"
        f"   commonf3x4 = 0\n"
    )
    proc = _run_paml([command, ctl.name], workdir, timeout)
    out = workdir / "yn"
    if not out.exists():
        raise RuntimeError(
            f"YN00 produced no output (exit {proc.returncode}): {proc.stdout.strip()[-300:]}"
        )
    text = out.read_text()
    # The output carries several methods; take the Yang & Nielsen (2000) block.
    marker = text.find(_YN00_SECTION)
    if marker == -1:
        raise RuntimeError("YN00 output contained no Yang & Nielsen section")
    match = _RE_YN00.search(text[marker:])
    if not match:
        raise RuntimeError("YN00 output contained no result line")
    return {key: _to_float(value) for key, value in match.groupdict().items()}


def _to_float(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


# --------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------


COLUMNS = [
    "orthogroup",
    "species",
    "gene_a",
    "gene_b",
    "n_seqs_alignment",
    "n_codons_alignment",
    "has_tree",
    "crosschecked",
    "n_codons_pair",
    "n_codons_a",
    "n_codons_b",
    "pct_id_a_in_b",
    "pct_id_b_in_a",
    "tree_dS",
    "tree_dN",
    "tree_omega",
    "pair_dS",
    "pair_dN",
    "pair_omega",
    "pair_t",
    "pair_S",
    "pair_N",
    "yn00_dS",
    "yn00_dN",
    "yn00_omega",
    "m0_omega",
    "m0_kappa",
    "m0_lnL",
    "dS_tree_over_pair",
]


def fmt(value) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return "NA"
    if isinstance(value, float):
        return f"{value:.6g}"
    return str(value)


def write_tsv(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({key: fmt(row.get(key)) for key in COLUMNS})


# --------------------------------------------------------------------------


def main() -> int:
    args = parse_args()
    msa_path = Path(args.msa)
    output_path = Path(args.output)
    orthogroup = args.orthogroup or msa_path.stem

    if args.icode not in PAML_ICODE_TO_NCBI:
        print(
            f"ERROR: --icode {args.icode} is not a PAML genetic code index "
            f"(expected one of {sorted(PAML_ICODE_TO_NCBI)}).",
            file=sys.stderr,
        )
        return 1
    codon_table = CodonTable.unambiguous_dna_by_id[PAML_ICODE_TO_NCBI[args.icode]]

    try:
        alignment = read_protein_alignment(msa_path, args.msa_format)
        cds_index = build_cds_index([Path(p) for p in args.cds])
        paralogs = read_manifest(Path(args.paralogs))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    protein_ids = [record.id for record in alignment]
    matched, missing = match_cds(protein_ids, cds_index)
    if missing:
        print(
            f"WARNING: {len(missing)} sequence(s) in {orthogroup} have no CDS match and "
            f"will be dropped (e.g. {missing[:3]}).",
            file=sys.stderr,
        )
        drop = set(missing)
        alignment = MultipleSeqAlignment([r for r in alignment if r.id not in drop])
        protein_ids = [pid for pid in protein_ids if pid not in drop]

    # A cheap guard before the codon alignment, which is the expensive step.
    if not select_reportable(paralogs, set(protein_ids)):
        total = sum(len(g) for g in paralogs.values())
        print(
            f"ERROR: No target species retains 2 or more paralogs in {orthogroup} after "
            f"CDS matching ({total} manifest gene(s) across {len(paralogs)} species). "
            f"Nothing to report.",
            file=sys.stderr,
        )
        return 1

    try:
        codon_alignment, unusable = build_codon_alignment_tolerantly(
            alignment, matched, codon_table
        )
    except Exception as exc:
        print(f"ERROR: Failed to build codon alignment for {orthogroup}: {exc}", file=sys.stderr)
        return 1

    # Everything the CDS files are needed for has happened; `matched` holds the
    # records themselves, so the open offset indexes can go.
    cds_index.close()

    if unusable:
        print(
            f"WARNING: {len(unusable)} sequence(s) in {orthogroup} do not reconcile with "
            f"their CDS and were dropped (e.g. {unusable[:3]}).",
            file=sys.stderr,
        )
        drop = set(unusable)
        alignment = MultipleSeqAlignment([r for r in alignment if r.id not in drop])
        protein_ids = [pid for pid in protein_ids if pid not in drop]

    # Re-derived after the drop: both the pair list and the gene tree have to
    # match the sequences that actually made it into the codon alignment.
    retained = set(protein_ids)
    reportable_by_species = select_reportable(paralogs, retained)
    if not reportable_by_species:
        print(
            f"ERROR: No target species retains 2 or more paralogs in {orthogroup} after "
            f"dropping {len(unusable)} sequence(s) that do not reconcile with their CDS.",
            file=sys.stderr,
        )
        return 1
    reportable = [gene for genes in reportable_by_species.values() for gene in genes]
    if len(protein_ids) < 2:
        print(f"ERROR: Fewer than 2 sequences remain in {orthogroup}.", file=sys.stderr)
        return 1

    id_map = make_short_id_map(protein_ids)
    reverse_id_map = {short: original for original, short in id_map.items()}
    codon_seqs = {record.id: str(record.seq) for record in codon_alignment}
    codon_rows = [(id_map[pid], codon_seqs[pid]) for pid in protein_ids if pid in codon_seqs]
    n_codons_alignment = len(codon_rows[0][1]) // 3 if codon_rows else 0

    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_root = Path(tempfile.mkdtemp(prefix=f"ks_{orthogroup}_"))
    tree_distances: dict[tuple[str, str], tuple[float, float]] = {}
    m0_stats: dict[str, float] = {}
    has_tree = False

    tree_path = Path(args.tree) if args.tree else None
    if tree_path and not args.skip_tree and tree_path.stat().st_size == 0:
        print(
            f"NOTE: No gene tree available for {orthogroup}; pairwise estimates only.",
            file=sys.stderr,
        )
        tree_path = None

    if tree_path and not args.skip_tree:
        try:
            tree = load_and_prepare_tree(tree_path, retained, id_map)
            tree_distances, m0_stats, edges, short_by_node = codeml_tree_estimates(
                codon_rows,
                tree,
                [id_map[gene] for gene in reportable],
                tmp_root / "tree",
                args.icode,
                args.codeml_method,
                args.codeml_command,
                args.codeml_timeout,
            )
            has_tree = True

            original_by_node = {
                node: reverse_id_map.get(short, short) for node, short in short_by_node.items()
            }
            for weight in ("dS", "dN"):
                newick = newick_from_branch_table(edges, original_by_node, weight)
                if newick is None:
                    print(
                        f"WARNING: Could not rebuild the {weight} tree for {orthogroup}.",
                        file=sys.stderr,
                    )
                    continue
                (output_path.parent / f"{orthogroup}_{weight}.nwk").write_text(newick + "\n")
        except Exception as exc:
            # A missing or unusable tree costs us the primary estimate but the
            # pairwise columns are still worth producing.
            print(
                f"WARNING: Tree-based estimate unavailable for {orthogroup} ({exc}). "
                f"Falling back to pairwise estimates only.",
                file=sys.stderr,
            )

    rows = []
    pairs = [
        (species, gene_a, gene_b)
        for species, genes in reportable_by_species.items()
        for gene_a, gene_b in itertools.combinations(genes, 2)
    ]

    crosscheck = set(pairs)
    if crosscheck_cap_applies(len(codon_rows), len(pairs), args.max_pairwise_pairs):
        # Seeded on the orthogroup so a rerun cross-checks the same pairs.
        crosscheck = set(
            random.Random(orthogroup).sample(pairs, args.max_pairwise_pairs)
        )
        covered = (
            "tree_dS still covers every pair"
            if has_tree
            else "the tree fit failed, so the remaining pairs get no estimate -- but "
            "running all of them would cost more than the fit they were checking"
        )
        print(
            f"NOTE: {orthogroup} has {len(pairs)} pairs; cross-checking a random "
            f"{args.max_pairwise_pairs} of them. {covered}.",
            file=sys.stderr,
        )

    # Phase one: every pair gets a row, with everything that costs nothing --
    # the tree-based dS included, since it was already computed above.
    pending: list[tuple[dict, str, str, str, str, str, str]] = []
    for species, gene_a, gene_b in pairs:
        short_a, short_b = id_map[gene_a], id_map[gene_b]
        aligned_a, aligned_b = codon_seqs[gene_a], codon_seqs[gene_b]
        pair_a, pair_b = strip_pair_gaps(aligned_a, aligned_b)
        id_a_in_b, id_b_in_a = pct_identity(aligned_a, aligned_b)
        wanted = (species, gene_a, gene_b) in crosscheck

        row: dict = {
            "orthogroup": orthogroup,
            "species": species,
            "gene_a": gene_a,
            "gene_b": gene_b,
            "n_seqs_alignment": len(codon_rows),
            "n_codons_alignment": n_codons_alignment,
            "has_tree": "yes" if has_tree else "no",
            "crosschecked": "yes" if wanted else "no",
            "n_codons_pair": len(pair_a) // 3,
            # Ungapped length of each sequence. n_codons_pair on its own cannot
            # separate a genuinely divergent pair from one where a sequence is a
            # fragment, and n_codons_alignment is the wrong denominator: it is
            # the union of every indel in the family, so a healthy orthogroup of
            # variable-length proteins looks as bad as a broken one. Coverage
            # against min(n_codons_a, n_codons_b) is the comparison that works.
            "n_codons_a": ungapped_codons(aligned_a),
            "n_codons_b": ungapped_codons(aligned_b),
            "pct_id_a_in_b": id_a_in_b * 100,
            "pct_id_b_in_a": id_b_in_a * 100,
            "m0_omega": m0_stats.get("omega"),
            "m0_kappa": m0_stats.get("kappa"),
            "m0_lnL": m0_stats.get("lnL"),
        }

        distance = tree_distances.get((short_a, short_b)) or tree_distances.get((short_b, short_a))
        if distance is not None:
            row["tree_dS"], row["tree_dN"] = distance
            if distance[0] and not math.isnan(distance[0]):
                row["tree_omega"] = distance[1] / distance[0]

        rows.append(row)
        if not wanted:
            continue
        if len(pair_a) < 3:
            print(
                f"WARNING: {gene_a} vs {gene_b} in {orthogroup} share no ungapped codons; "
                f"pairwise estimates skipped.",
                file=sys.stderr,
            )
            row["crosschecked"] = "no"
            continue
        pending.append((row, short_a, pair_a, short_b, pair_b, gene_a, gene_b))

    # Phase two: the cross-checks, which are two external processes per pair and
    # the only expensive thing left. Threads rather than processes because the
    # work happens in CODEML and YN00, not in this interpreter; each pair writes
    # into its own directory, so they do not interact.
    def crosscheck_pair(job) -> None:
        row, short_a, pair_a, short_b, pair_b, gene_a, gene_b = job
        pair_dir = tmp_root / f"pair_{short_a}_{short_b}"
        if not args.skip_pairwise:
            try:
                result = codeml_pairwise_estimate(
                    short_a, pair_a, short_b, pair_b, pair_dir / "codeml",
                    args.icode, args.codeml_command, args.codeml_timeout,
                )
                row["pair_dS"] = result["dS"]
                row["pair_dN"] = result["dN"]
                row["pair_omega"] = result["omega"]
                row["pair_t"] = result["t"]
                row["pair_S"] = result["S"]
                row["pair_N"] = result["N"]
            except Exception as exc:
                print(
                    f"WARNING: pairwise CODEML failed for {gene_a} vs {gene_b} "
                    f"in {orthogroup}: {exc}",
                    file=sys.stderr,
                )
        if not args.skip_yn00:
            try:
                result = yn00_pairwise_estimate(
                    short_a, pair_a, short_b, pair_b, pair_dir / "yn00",
                    args.icode, args.yn00_command, args.codeml_timeout,
                )
                row["yn00_dS"] = result["dS"]
                row["yn00_dN"] = result["dN"]
                row["yn00_omega"] = result["omega"]
            except Exception as exc:
                print(
                    f"WARNING: YN00 failed for {gene_a} vs {gene_b} in {orthogroup}: {exc}",
                    file=sys.stderr,
                )

    if pending:
        if args.threads > 1:
            with ThreadPoolExecutor(max_workers=args.threads) as pool:
                list(pool.map(crosscheck_pair, pending))
        else:
            for job in pending:
                crosscheck_pair(job)

    for row in rows:
        tree_ds, pair_ds = row.get("tree_dS"), row.get("pair_dS")
        if (
            isinstance(tree_ds, float)
            and isinstance(pair_ds, float)
            and not math.isnan(tree_ds)
            and not math.isnan(pair_ds)
            and pair_ds > 0
        ):
            row["dS_tree_over_pair"] = tree_ds / pair_ds

    try:
        write_tsv(rows, output_path)
    except Exception as exc:
        print(f"ERROR: Failed to write {output_path}: {exc}", file=sys.stderr)
        return 1
    finally:
        if args.keep_temp:
            print(f"NOTE: CODEML working files kept in {tmp_root}", file=sys.stderr)
        else:
            shutil.rmtree(tmp_root, ignore_errors=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
