#!/usr/bin/env python3
"""Reduce a proteome to one representative isoform per gene.

Two schemes, because the annotation sources differ in where the gene identifier
lives:

  default   The protein defline carries it (Ensembl-style `gene:ENSDARG...`).
            Selection reads the protein FASTA and keeps the longest protein per
            gene. `--regex` sets the pattern.

  --refseq  RefSeq protein deflines carry no gene identifier at all:

                >NP_001001398.2 fibroblast growth factor 6a precursor [Danio rerio]

            but the CDS FASTA carries both the gene and the protein it encodes:

                >lcl|NC_133200.1_cds_NP_001001398.2_93141 [gene=fgf6a]
                 [db_xref=GeneID:30560] [protein_id=NP_001001398.2] ...

            so selection is driven from the CDS file: group by GeneID, keep the
            longest CDS, and pull its protein out of the .faa by protein_id.

`--all` applies either scheme across a directory of per-species subdirectories.

When CDS is supplied the outputs are a matched pair whose record IDs are
identical -- the bare protein accession, no description -- which is the form the
pipeline requires, since `ks.py` matches CDS to protein on exact ID equality.
CDS is written verbatim, trailing stop codon included.

Examples:
    longest_isoform.py --input proteome.faa                       # Ensembl-style
    longest_isoform.py --refseq --input x_protein.faa --cds x_cds_from_genomic.fna
    longest_isoform.py --refseq --all ./fish_data
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# NCBI writes CDS metadata as bracketed key=value fields. A value never contains
# ']', so this splits them safely even when one holds commas or dots.
_FIELD_RE = re.compile(r"\[([^\[\]=]+)=([^\]]*)\]")
# db_xref packs several sources into one field: [db_xref=GeneID:30560,ZFIN:...]
_GENEID_RE = re.compile(r"\bGeneID:(\d+)\b")
# Fallback for records missing [protein_id=...]:
# lcl|NC_133200.1_cds_NP_001001398.2_93141 -> NP_001001398.2
_LCL_PROTEIN_RE = re.compile(r"_cds_(.+)_\d+$")

DEFAULT_REGEX = r"gene:\S+"
DEFAULT_PROTEIN_GLOB = "*_protein.faa"
DEFAULT_CDS_GLOB = "*_cds_from_genomic.fna"


# --------------------------------------------------------------------------- #
# Defline parsing
# --------------------------------------------------------------------------- #

def parse_bracket_fields(description: str) -> dict[str, str]:
    """Bracketed key=value fields from an NCBI CDS defline.

    Repeated db_xref fields accumulate; everything else keeps its first value.
    """
    fields: dict[str, str] = {}
    for key, value in _FIELD_RE.findall(description):
        if key == "db_xref":
            fields["db_xref"] = f"{fields.get('db_xref', '')},{value}"
        else:
            fields.setdefault(key, value)
    return fields


def parse_refseq_cds_defline(description: str, record_id: str) -> tuple[str | None, str | None, str | None]:
    """Return (gene_key, protein_id, gene_symbol) for one RefSeq CDS record.

    gene_key prefers the GeneID, which is stable and unambiguous; a record
    without one falls back to its locus tag and then its gene symbol rather than
    being dropped, since either still groups isoforms correctly within one
    assembly. protein_id is absent for pseudogenes and non-coding features.
    """
    fields = parse_bracket_fields(description)

    protein_id = fields.get("protein_id")
    if not protein_id:
        match = _LCL_PROTEIN_RE.search(record_id)
        protein_id = match.group(1) if match else None

    symbol = fields.get("gene")
    geneid = _GENEID_RE.search(fields.get("db_xref", ""))
    if geneid:
        gene_key = f"GeneID:{geneid.group(1)}"
    elif fields.get("locus_tag"):
        gene_key = f"locus_tag:{fields['locus_tag']}"
    else:
        gene_key = f"gene:{symbol}" if symbol else None

    return gene_key, protein_id, symbol


def gene_key_from_regex(description: str, pattern: re.Pattern) -> str | None:
    """First regex match in a protein defline, or None."""
    match = pattern.search(description)
    return match.group(0) if match else None


# --------------------------------------------------------------------------- #
# Selection
# --------------------------------------------------------------------------- #

@dataclass
class Isoform:
    """The best isoform seen so far for one gene."""

    gene_key: str
    protein_id: str = ""
    symbol: str | None = None
    protein: SeqRecord | None = None
    cds: Seq | None = None
    n_isoforms: int = 0

    def length(self) -> int:
        """The quantity being maximised: CDS length when known, else protein."""
        if self.cds is not None:
            return len(self.cds)
        return len(self.protein.seq) if self.protein else 0


def load_fasta_index(path: Path) -> dict[str, SeqRecord]:
    """Index a FASTA by bare record id, warning on (and keeping) the first duplicate."""
    index: dict[str, SeqRecord] = {}
    with open(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in index:
                print(f"WARNING: duplicate ID '{record.id}' in {path}; keeping the first.", file=sys.stderr)
                continue
            index[record.id] = record
    if not index:
        raise ValueError(f"No sequences read from: {path}")
    return index


def select_refseq(cds_paths: list[Path], proteins: dict[str, SeqRecord], stats: dict[str, int]) -> dict[str, Isoform]:
    """Longest CDS per gene, restricted to isoforms whose protein exists.

    Filtering on protein availability during selection (rather than after) means
    a gene is not lost when its longest isoform happens to have no protein
    record -- the next-longest one that resolves wins instead.
    """
    best: dict[str, Isoform] = {}
    for path in cds_paths:
        with open(path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                stats["records"] += 1
                gene_key, protein_id, symbol = parse_refseq_cds_defline(record.description, record.id)
                if gene_key is None:
                    stats["no_gene_key"] += 1
                    continue
                if protein_id is None:
                    stats["no_protein_id"] += 1
                    continue

                current = best.setdefault(gene_key, Isoform(gene_key, symbol=symbol))
                current.n_isoforms += 1
                if protein_id not in proteins:
                    stats["protein_missing"] += 1
                    continue
                if current.cds is None or len(record.seq) > len(current.cds):
                    current.protein_id = protein_id
                    current.protein = proteins[protein_id]
                    current.cds = record.seq
                    current.symbol = symbol or current.symbol

    return {key: iso for key, iso in best.items() if iso.cds is not None}


def select_by_regex(
    protein_path: Path, pattern: re.Pattern, cds: dict[str, SeqRecord] | None, stats: dict[str, int]
) -> dict[str, Isoform]:
    """Longest protein per gene, with the gene id taken from the protein defline."""
    best: dict[str, Isoform] = {}
    with open(protein_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            stats["records"] += 1
            gene_key = gene_key_from_regex(record.description, pattern)
            if gene_key is None:
                stats["no_gene_key"] += 1
                continue

            current = best.setdefault(gene_key, Isoform(gene_key))
            current.n_isoforms += 1
            if cds is not None and record.id not in cds:
                stats["cds_missing"] += 1
                continue
            if current.protein is None or len(record.seq) > len(current.protein.seq):
                current.protein_id = record.id
                current.protein = record
                current.cds = cds[record.id].seq if cds is not None else None

    return {key: iso for key, iso in best.items() if iso.protein is not None}


def drop_duplicate_accessions(best: dict[str, Isoform], stats: dict[str, int]) -> dict[str, Isoform]:
    """Keep one gene per protein accession, preferring the longest.

    The same accession can legitimately appear under two gene keys on alternate
    loci or PAR regions, but `ks.py` treats a duplicate ID as a hard error.
    """
    seen: dict[str, str] = {}
    for key in sorted(best, key=lambda k: (-best[k].length(), k)):
        pid = best[key].protein_id
        if pid in seen:
            stats["duplicate_protein"] += 1
            del best[key]
        else:
            seen[pid] = key
    return best


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #

def add_suffix(path: Path, suffix: str) -> Path:
    """`x/GCF_1.2_hap1.1_protein.faa` -> `x/GCF_1.2_hap1.1_protein_longest.faa`.

    Only the final extension is separated, so the dots inside RefSeq accessions
    survive untouched.
    """
    return path.with_name(f"{path.stem}_{suffix}{path.suffix}")


def translates_to_protein(iso: Isoform, table: int) -> bool:
    """True when the chosen CDS translates to exactly its protein."""
    if iso.cds is None or iso.protein is None:
        return True
    try:
        observed = str(iso.cds.translate(table=table)).rstrip("*")
    except Exception:
        return False
    return observed == str(iso.protein.seq).rstrip("*")


def write_outputs(
    best: dict[str, Isoform], out_protein: Path, out_cds: Path | None,
    report: Path | None, table: int, drop_bad: bool,
) -> tuple[int, int]:
    """Write the protein FASTA, optional CDS FASTA and optional report.

    Returns (pairs written, translation mismatches seen).
    """
    proteins: list[SeqRecord] = []
    cds_records: list[SeqRecord] = []
    rows: list[tuple] = []
    mismatches = 0

    for key, iso in sorted(best.items()):
        ok = translates_to_protein(iso, table)
        if not ok:
            mismatches += 1
            if drop_bad:
                continue
        proteins.append(SeqRecord(seq=iso.protein.seq, id=iso.protein_id, description=""))
        if iso.cds is not None:
            cds_records.append(SeqRecord(seq=iso.cds, id=iso.protein_id, description=""))
        rows.append((key, iso.symbol or "NA", iso.protein_id, iso.length(), iso.n_isoforms,
                     "yes" if ok else "no"))

    out_protein.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(proteins, out_protein, "fasta")
    if out_cds is not None and cds_records:
        out_cds.parent.mkdir(parents=True, exist_ok=True)
        SeqIO.write(cds_records, out_cds, "fasta")

    if report is not None:
        report.parent.mkdir(parents=True, exist_ok=True)
        with open(report, "w", newline="") as handle:
            handle.write("gene_key\tgene_symbol\tprotein_id\tlength\tn_isoforms\ttranslation_ok\n")
            for row in rows:
                handle.write("\t".join(str(field) for field in row) + "\n")

    return len(proteins), mismatches


# --------------------------------------------------------------------------- #
# Drivers
# --------------------------------------------------------------------------- #

def process_one(
    protein_path: Path, cds_paths: list[Path], out_protein: Path, out_cds: Path | None,
    report: Path | None, args: argparse.Namespace,
) -> int:
    """Select isoforms for one species. Returns the number of pairs written."""
    stats = dict.fromkeys(
        ("records", "no_gene_key", "no_protein_id", "protein_missing", "cds_missing", "duplicate_protein"), 0
    )

    if args.refseq:
        proteins = load_fasta_index(protein_path)
        best = select_refseq(cds_paths, proteins, stats)
    else:
        cds_index = None
        if cds_paths:
            cds_index = {}
            for path in cds_paths:
                cds_index.update(load_fasta_index(path))
        best = select_by_regex(protein_path, re.compile(args.regex), cds_index, stats)

    if not best:
        raise ValueError(f"No gene could be resolved from {protein_path}")

    best = drop_duplicate_accessions(best, stats)
    written, mismatches = write_outputs(
        best, out_protein, out_cds, report, args.table, args.drop_untranslatable
    )

    print(f"Read {stats['records']:,} record(s); resolved {len(best):,} gene(s); "
          f"wrote {written:,}.", file=sys.stderr)
    for label, key in (
        ("no gene identifier", "no_gene_key"),
        ("no protein_id (pseudogene/non-coding)", "no_protein_id"),
        ("protein absent from the .faa", "protein_missing"),
        ("no matching CDS record", "cds_missing"),
        ("dropped as duplicate accessions", "duplicate_protein"),
    ):
        if stats[key]:
            print(f"  {stats[key]:,}: {label}.", file=sys.stderr)
    if mismatches:
        verb = "dropped" if args.drop_untranslatable else "kept"
        print(f"  WARNING: {mismatches:,} CDS did not translate to its protein "
              f"under table {args.table} ({verb}).", file=sys.stderr)
        if args.require_translation:
            raise ValueError("--require-translation was set")
    return written


def find_one(directory: Path, pattern: str, what: str) -> Path | None:
    """The single file matching `pattern`, or None with a reason on stderr.

    Zero or several matches is a download problem; guessing would silently pair
    the wrong files.
    """
    hits = sorted(directory.glob(pattern))
    if len(hits) != 1:
        print(f"SKIP {directory.name}: {len(hits)} {what} match '{pattern}'; need exactly 1.",
              file=sys.stderr)
        return None
    return hits[0]


def process_directory(root: Path, args: argparse.Namespace) -> int:
    """Apply the selected scheme to every subdirectory of `root`."""
    done = skipped = failed = 0
    for directory in sorted(p for p in root.iterdir() if p.is_dir()):
        protein = find_one(directory, args.protein_glob, "protein file(s)")
        cds = find_one(directory, args.cds_glob, "CDS file(s)") if args.refseq or args.cds_glob else None
        if protein is None or (args.refseq and cds is None):
            skipped += 1
            continue

        print(f"=== {directory.name} ===", file=sys.stderr)
        try:
            # One malformed download must not abandon the other species.
            process_one(
                protein, [cds] if cds else [],
                add_suffix(protein, args.suffix),
                add_suffix(cds, args.suffix) if cds else None,
                directory / f"{directory.name}_isoforms.tsv",
                args,
            )
            done += 1
        except (OSError, ValueError) as exc:
            print(f"FAILED {directory.name}: {exc}", file=sys.stderr)
            failed += 1

    print(f"\nProcessed {done} species ({skipped} skipped, {failed} failed).", file=sys.stderr)
    return 1 if failed else 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Reduce a proteome to one representative isoform per gene.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--input", type=Path, help="Protein FASTA for a single species.")
    source.add_argument("--all", type=Path, metavar="ROOT",
                        help="Directory with one subdirectory per species.")

    parser.add_argument("--refseq", action="store_true",
                        help="Drive selection from the CDS file using RefSeq's "
                             "[db_xref=GeneID:...] fields instead of the protein defline.")
    parser.add_argument("--regex", default=DEFAULT_REGEX,
                        help=f"Gene-id pattern in the protein defline, ignored with "
                             f"--refseq (default: {DEFAULT_REGEX!r}).")
    parser.add_argument("--cds", nargs="+", type=Path, default=[],
                        help="CDS FASTA(s); required with --refseq, optional otherwise.")

    parser.add_argument("--out-proteins", type=Path, help="Output protein FASTA.")
    parser.add_argument("--out-cds", type=Path, help="Output CDS FASTA.")
    parser.add_argument("--report", type=Path, help="Optional TSV of the isoform chosen per gene.")
    parser.add_argument("--suffix", default="longest",
                        help="Inserted before the final extension when an output "
                             "path is not given (default: longest).")

    parser.add_argument("--protein-glob", default=DEFAULT_PROTEIN_GLOB,
                        help=f"--all: protein filename pattern (default: {DEFAULT_PROTEIN_GLOB}).")
    parser.add_argument("--cds-glob", default=DEFAULT_CDS_GLOB,
                        help=f"--all: CDS filename pattern (default: {DEFAULT_CDS_GLOB}).")

    parser.add_argument("--table", type=int, default=1,
                        help="NCBI genetic code used to verify translations (default: 1).")
    parser.add_argument("--require-translation", action="store_true",
                        help="Exit non-zero if any chosen CDS does not translate to its protein.")
    parser.add_argument("--drop-untranslatable", action="store_true",
                        help="Omit genes whose chosen CDS does not translate to its protein.")

    args = parser.parse_args()
    if args.refseq and not (args.all or args.cds):
        parser.error("--refseq needs --cds (or --all, which discovers it)")
    return args


def main() -> int:
    args = parse_args()
    try:
        if args.all:
            return process_directory(args.all, args)
        process_one(
            args.input, args.cds,
            args.out_proteins or add_suffix(args.input, args.suffix),
            args.out_cds or (add_suffix(args.cds[0], args.suffix) if args.cds else None),
            args.report, args,
        )
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
