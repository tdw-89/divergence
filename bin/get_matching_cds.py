#!/usr/bin/env python3
"""
Find CDS nucleotide sequences that correspond to proteins in a protein FASTA file.

A CDS record is accepted only when BOTH conditions hold:
  1. The protein ID appears somewhere in the CDS defline (any field, not just 'gene:').
     Matching is done on the base ID (version suffix stripped) for robustness.
  2. Translating the CDS nucleotide sequence yields exactly the protein's amino
     acid sequence (standard genetic code, trailing stop codon ignored).

Output is a FASTA file with the matching CDS nucleotide sequences and only the
original protein ID as the defline — no additional description.

Usage:
    python get_matching_cds.py <protein_fasta> <cds_fasta> <output_fasta>
"""

import argparse
import re
import sys

import Bio.SeqIO
from Bio.SeqRecord import SeqRecord

# Matches any alphanumeric+dot token (captures IDs embedded after colons, pipes, etc.)
_TOKEN_RE = re.compile(r'[\w.]+')


def load_proteins(protein_fasta: str) -> dict[str, tuple[str, str]]:
    """
    Parse protein FASTA.
    Returns: {base_id: (full_versioned_id, aa_sequence_without_trailing_stop)}
    """
    proteins: dict[str, tuple[str, str]] = {}
    for record in Bio.SeqIO.parse(protein_fasta, "fasta"):
        full_id = record.id
        base_id = full_id.split(".")[0]
        aa_seq = str(record.seq).rstrip("*")
        if base_id in proteins:
            print(
                f"Warning: duplicate protein base ID '{base_id}'; keeping first.",
                file=sys.stderr,
            )
        else:
            proteins[base_id] = (full_id, aa_seq)
    return proteins

def candidates_from_defline(defline: str, proteins: dict) -> list[str]:
    """
    Return base protein IDs that appear anywhere in the CDS defline.

    Splits the defline into alphanumeric tokens (so 'gene:ENSDARG00000086269.4'
    yields the token 'ENSDARG00000086269.4'), strips version suffixes, and
    checks each against the protein ID lookup table.
    """
    seen: set[str] = set()
    for token in _TOKEN_RE.findall(defline):
        base = token.split(".")[0]
        if base in proteins and base not in seen:
            seen.add(base)
    return list(seen)

def translate_cds(seq) -> str:
    """Translate a CDS sequence (standard genetic code); strip any trailing stop codon."""
    return str(seq.translate(table=1)).rstrip("*")

def select_best(candidates: list[tuple[SeqRecord, str]]) -> SeqRecord:
    """
    Pick the best CDS record from a list of (SeqRecord, original_defline) tuples.

    Priority:
      1. Longest sequence.
      2. On length tie, prefer records whose defline contains 'protein_coding'.
      3. On further tie, keep the first encountered (stable sort).
    """
    def sort_key(item: tuple[SeqRecord, str]) -> tuple[int, int]:
        record, defline = item
        return (-len(str(record.seq)), 0 if "protein_coding" in defline else 1)

    best_record, _ = sorted(candidates, key=sort_key)[0]
    return best_record


def find_matching_cds(protein_fasta: str, cds_fasta: str, output_fasta: str) -> None:
    proteins = load_proteins(protein_fasta)
    if not proteins:
        print("Error: no protein records found in the protein FASTA.", file=sys.stderr)
        sys.exit(1)

    print(f"Loaded {len(proteins)} protein records.", file=sys.stderr)

    # full_id -> list of (SeqRecord with protein id, original defline)
    matches: dict[str, list[tuple[SeqRecord, str]]] = {}
    n_candidates = 0
    n_seq_mismatches = 0

    for record in Bio.SeqIO.parse(cds_fasta, "fasta"):
        matched_ids = candidates_from_defline(record.description, proteins)
        if not matched_ids:
            continue

        n_candidates += 1
        translated = translate_cds(record.seq)

        for base_id in matched_ids:
            full_id, aa_seq = proteins[base_id]
            if translated == aa_seq:
                entry = (SeqRecord(seq=record.seq, id=full_id, description=""), record.description)
                matches.setdefault(full_id, []).append(entry)
            else:
                n_seq_mismatches += 1

    print(
        f"Checked {n_candidates} CDS candidates (defline match); "
        f"{n_seq_mismatches} discarded after sequence verification.",
        file=sys.stderr,
    )

    unmatched = [full_id for base_id, (full_id, _) in proteins.items() if full_id not in matches]
    if unmatched:
        print(
            f"Warning: {len(unmatched)} protein(s) had no matching CDS sequence:",
            file=sys.stderr,
        )
        for pid in unmatched:
            print(f"  {pid}", file=sys.stderr)

    output_records = [select_best(candidates) for candidates in matches.values()]

    Bio.SeqIO.write(output_records, output_fasta, "fasta")
    print(
        f"Wrote {len(output_records)} CDS records (one per unique protein ID) to '{output_fasta}'.",
        file=sys.stderr,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Find CDS nucleotide sequences that match a given protein by both "
            "ID presence in the CDS defline and exact translated amino acid sequence."
        )
    )
    parser.add_argument("protein_fasta", help="Input protein FASTA file")
    parser.add_argument("cds_fasta", help="Input CDS nucleotide FASTA file")
    parser.add_argument("output_fasta", help="Output FASTA file for matching CDS sequences")
    args = parser.parse_args()

    find_matching_cds(args.protein_fasta, args.cds_fasta, args.output_fasta)


if __name__ == "__main__":
    main()

