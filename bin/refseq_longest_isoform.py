#!/usr/bin/env python3
"""Pick one representative isoform per gene from an NCBI RefSeq download.

RefSeq splits the information this needs across two files. The protein FASTA
carries no gene identifier at all:

    >NP_001001398.2 fibroblast growth factor 6a precursor [Danio rerio]

while the CDS FASTA carries both the gene and the protein it encodes:

    >lcl|NC_133200.1_cds_NP_001001398.2_93141 [gene=fgf6a]
     [db_xref=GeneID:30560] [protein=...] [protein_id=NP_001001398.2] ...

So isoform selection has to be driven from the CDS file: group CDS records by
GeneID, take the longest, and pull the matching protein out of the .faa by its
protein_id. A regex over the protein deflines (what `longest_isoform.py` does
for Ensembl-style headers) cannot work here.

Writes a protein FASTA and a CDS FASTA whose record IDs are identical -- the
bare versioned protein accession, no description -- which is the form the
pipeline requires, since `ks.py` matches CDS to protein on exact ID equality.
The CDS is written exactly as NCBI provides it, trailing stop codon included;
`Bio.codonalign` accepts that and the rest of the pipeline assumes it.

Usage:
    refseq_longest_isoform.py --cds GCF_x_cds_from_genomic.fna \\
        --proteins GCF_x_protein.faa \\
        --out-proteins danio.faa --out-cds danio.fna
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# NCBI writes CDS metadata as bracketed key=value fields on the defline. A value
# never contains ']', so this is a safe way to pull them apart.
_FIELD_RE = re.compile(r"\[([^\[\]=]+)=([^\]]*)\]")
# db_xref packs several sources into one field: [db_xref=GeneID:30560,ZFIN:ZDB-...]
_GENEID_RE = re.compile(r"\bGeneID:(\d+)\b")
# Fallback for records missing [protein_id=...]:
# lcl|NC_133200.1_cds_NP_001001398.2_93141 -> NP_001001398.2
_LCL_PROTEIN_RE = re.compile(r"_cds_(.+)_\d+$")


def parse_cds_defline(description: str, record_id: str) -> tuple[str | None, str | None, str | None]:
    """Return (gene_key, protein_id, gene_symbol) for one CDS defline.

    `gene_key` prefers the GeneID, which is stable and unambiguous. When a
    record carries no GeneID we fall back to the locus tag and then the gene
    symbol rather than dropping the record, since either still groups isoforms
    of the same gene correctly within a single assembly.
    """
    fields: dict[str, str] = {}
    for key, value in _FIELD_RE.findall(description):
        # Keep the first occurrence; NCBI does not repeat these except db_xref,
        # which we accumulate below.
        if key == "db_xref":
            fields["db_xref"] = fields.get("db_xref", "") + "," + value
        else:
            fields.setdefault(key, value)

    protein_id = fields.get("protein_id")
    if not protein_id:
        match = _LCL_PROTEIN_RE.search(record_id)
        protein_id = match.group(1) if match else None

    gene_symbol = fields.get("gene")

    gene_key = None
    geneid = _GENEID_RE.search(fields.get("db_xref", ""))
    if geneid:
        gene_key = f"GeneID:{geneid.group(1)}"
    elif fields.get("locus_tag"):
        gene_key = f"locus_tag:{fields['locus_tag']}"
    elif gene_symbol:
        gene_key = f"gene:{gene_symbol}"

    return gene_key, protein_id, gene_symbol


def load_proteins(path: Path) -> dict[str, SeqRecord]:
    """Index the protein FASTA by its bare record ID (the accession, versioned)."""
    proteins: dict[str, SeqRecord] = {}
    with open(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            if record.id in proteins:
                print(
                    f"WARNING: duplicate protein ID '{record.id}' in {path}; keeping the first.",
                    file=sys.stderr,
                )
                continue
            proteins[record.id] = record
    if not proteins:
        raise ValueError(f"No sequences read from protein FASTA: {path}")
    return proteins


class Selection:
    """The best CDS seen so far for one gene."""

    __slots__ = ("protein_id", "gene_symbol", "seq", "n_isoforms", "n_usable")

    def __init__(self, protein_id: str, gene_symbol: str | None) -> None:
        self.protein_id = protein_id
        self.gene_symbol = gene_symbol
        self.seq: Seq | None = None
        self.n_isoforms = 0
        self.n_usable = 0


def select_longest(
    cds_paths: list[Path], proteins: dict[str, SeqRecord]
) -> tuple[dict[str, Selection], dict[str, int]]:
    """Stream the CDS file(s), keeping the longest CDS per gene.

    Only candidates whose protein_id is present in the .faa are eligible, so a
    gene is not lost when its longest isoform happens to have no protein record
    -- the next-longest one that does resolve wins instead.
    """
    best: dict[str, Selection] = {}
    stats = {
        "records": 0,
        "no_gene_key": 0,
        "no_protein_id": 0,
        "protein_missing": 0,
        "duplicate_protein": 0,
    }

    for path in cds_paths:
        with open(path) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                stats["records"] += 1
                gene_key, protein_id, gene_symbol = parse_cds_defline(
                    record.description, record.id
                )
                if gene_key is None:
                    stats["no_gene_key"] += 1
                    continue
                if protein_id is None:
                    # Pseudogenes and non-coding features have no protein.
                    stats["no_protein_id"] += 1
                    continue

                current = best.get(gene_key)
                if current is None:
                    current = Selection(protein_id, gene_symbol)
                    best[gene_key] = current
                current.n_isoforms += 1

                if protein_id not in proteins:
                    stats["protein_missing"] += 1
                    continue
                current.n_usable += 1

                if current.seq is None or len(record.seq) > len(current.seq):
                    current.protein_id = protein_id
                    current.gene_symbol = gene_symbol or current.gene_symbol
                    current.seq = record.seq

    # Genes where no isoform resolved to a protein are not representable.
    unusable = [key for key, sel in best.items() if sel.seq is None]
    for key in unusable:
        del best[key]

    # The same protein accession can legitimately appear under two gene keys on
    # alternate loci or PAR regions; emitting it twice would break the pipeline's
    # requirement that IDs be unique, so keep the longest and drop the rest.
    by_protein: dict[str, str] = {}
    for key in sorted(best, key=lambda k: (-len(best[k].seq), k)):
        pid = best[key].protein_id
        if pid in by_protein:
            stats["duplicate_protein"] += 1
            del best[key]
        else:
            by_protein[pid] = key

    return best, stats


def validate_translation(selection: Selection, protein: SeqRecord, table: int) -> bool:
    """True when the CDS translates to exactly the protein sequence."""
    expected = str(protein.seq).rstrip("*")
    try:
        observed = str(selection.seq.translate(table=table)).rstrip("*")
    except Exception:
        return False
    return observed == expected


def main() -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Select the longest CDS isoform per gene from an NCBI RefSeq CDS FASTA "
            "and emit matched protein and CDS files with identical record IDs."
        )
    )
    parser.add_argument(
        "--cds", required=True, nargs="+", type=Path,
        help="RefSeq *_cds_from_genomic.fna file(s) for one species.",
    )
    parser.add_argument(
        "--proteins", required=True, type=Path,
        help="RefSeq *_protein.faa file for the same assembly.",
    )
    parser.add_argument("--out-proteins", required=True, type=Path, help="Output protein FASTA.")
    parser.add_argument("--out-cds", required=True, type=Path, help="Output CDS FASTA.")
    parser.add_argument(
        "--report", type=Path,
        help="Optional TSV listing the isoform chosen for each gene.",
    )
    parser.add_argument(
        "--table", type=int, default=1,
        help="NCBI genetic code table used to verify translations (default: 1).",
    )
    parser.add_argument(
        "--require-translation", action="store_true",
        help="Exit non-zero if any chosen CDS does not translate to its protein. "
        "Without this, mismatches are reported but the pair is still written.",
    )
    parser.add_argument(
        "--drop-untranslatable", action="store_true",
        help="Omit genes whose chosen CDS does not translate to its protein.",
    )
    args = parser.parse_args()

    try:
        proteins = load_proteins(args.proteins)
        best, stats = select_longest(args.cds, proteins)
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if not best:
        print("ERROR: no gene could be resolved to a protein/CDS pair.", file=sys.stderr)
        return 1

    protein_records: list[SeqRecord] = []
    cds_records: list[SeqRecord] = []
    report_rows: list[tuple[str, str, str, int, int, str]] = []
    n_mismatch = 0

    for gene_key, sel in sorted(best.items()):
        protein = proteins[sel.protein_id]
        ok = validate_translation(sel, protein, args.table)
        if not ok:
            n_mismatch += 1
            if args.drop_untranslatable:
                continue

        protein_records.append(
            SeqRecord(seq=protein.seq, id=sel.protein_id, description="")
        )
        cds_records.append(SeqRecord(seq=sel.seq, id=sel.protein_id, description=""))
        report_rows.append(
            (
                gene_key,
                sel.gene_symbol or "NA",
                sel.protein_id,
                len(sel.seq),
                sel.n_isoforms,
                "yes" if ok else "no",
            )
        )

    args.out_proteins.parent.mkdir(parents=True, exist_ok=True)
    args.out_cds.parent.mkdir(parents=True, exist_ok=True)
    SeqIO.write(protein_records, args.out_proteins, "fasta")
    SeqIO.write(cds_records, args.out_cds, "fasta")

    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        with open(args.report, "w", newline="") as handle:
            handle.write("gene_key\tgene_symbol\tprotein_id\tcds_length\tn_isoforms\ttranslation_ok\n")
            for row in report_rows:
                handle.write("\t".join(str(field) for field in row) + "\n")

    print(
        f"Read {stats['records']:,} CDS record(s); resolved {len(best):,} gene(s); "
        f"wrote {len(protein_records):,} protein/CDS pair(s).",
        file=sys.stderr,
    )
    for label, key in (
        ("no gene identifier", "no_gene_key"),
        ("no protein_id (pseudogene/non-coding)", "no_protein_id"),
        ("protein_id absent from the .faa", "protein_missing"),
        ("genes dropped as duplicate protein accessions", "duplicate_protein"),
    ):
        if stats[key]:
            print(f"  {stats[key]:,} CDS record(s): {label}.", file=sys.stderr)
    if n_mismatch:
        verb = "dropped" if args.drop_untranslatable else "kept"
        print(
            f"  WARNING: {n_mismatch:,} chosen CDS did not translate to its protein "
            f"under table {args.table} ({verb}).",
            file=sys.stderr,
        )
        if args.require_translation:
            print("ERROR: --require-translation was set.", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
