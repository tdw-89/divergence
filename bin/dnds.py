#!/usr/bin/env python3

"""Run YN00 dN/dS estimation from protein MSA + nucleotide FASTA.

Step 1: map MSA protein IDs to nucleotide records from a FASTA reference.
Step 2: build codon alignment using Bio.codonalign.build.
Step 3: run Bio.Phylo.PAML.yn00 and write pairwise results to TSV.
"""

from __future__ import annotations

import argparse
import csv
import math
import shutil
import sys
import tempfile
from pathlib import Path

from Bio import AlignIO, SeqIO
from Bio.codonalign import build
from Bio.Phylo.PAML import yn00
from Bio.SeqRecord import SeqRecord


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Run pairwise YN00 dN/dS from protein MSA and nucleotide FASTA."
	)
	parser.add_argument(
		"-m",
		"--msa",
		required=True,
		help="Path to the protein MSA file.",
	)
	parser.add_argument(
		"-f",
		"--fasta",
		required=True,
		help="Path to the nucleotide FASTA reference file.",
	)
	parser.add_argument(
		"--msa-format",
		default="fasta",
		help="Biopython alignment format for the protein MSA (default: fasta).",
	)
	parser.add_argument(
		"--yn00-command",
		default="yn00",
		help="Command used to run the external YN00 executable (default: yn00).",
	)
	parser.add_argument(
		"-o",
		"--output",
		required=True,
		help="Output TSV path for pairwise dN/dS results.",
	)
	parser.add_argument(
		"--yn00-verbose",
		action="store_true",
		help="Run yn00 in verbose mode (useful for debugging failures).",
	)
	return parser.parse_args()


def read_protein_alignment(msa_path: Path, msa_format: str):
	alignment = AlignIO.read(str(msa_path), msa_format)
	ids = [record.id for record in alignment]

	seen: set[str] = set()
	duplicate_ids: list[str] = []
	for seq_id in ids:
		if seq_id in seen:
			duplicate_ids.append(seq_id)
		else:
			seen.add(seq_id)
	if duplicate_ids:
		duplicates_str = ", ".join(sorted(set(duplicate_ids)))
		raise ValueError(
			f"Protein MSA has duplicate sequence IDs: {duplicates_str}. "
			"IDs must be unique to unambiguously map nucleotide sequences."
		)
	return alignment, ids


def build_nucleotide_index(fasta: Path) -> dict[str, SeqRecord]:
	index: dict[str, SeqRecord] = {}
	for record in SeqIO.parse(str(fasta), "fasta"):
		if record.id in index:
			raise ValueError(
				f"Duplicate nucleotide FASTA ID found: {record.id}. "
				"IDs must be unique for mapping."
			)
		index[record.id] = record
	return index


def find_matching_nucleotides(
	protein_ids: list[str],
	nucleotide_index: dict[str, SeqRecord],
) -> tuple[list[SeqRecord], list[str]]:
	matched: list[SeqRecord] = []
	missing: list[str] = []

	for protein_id in protein_ids:
		nucleotide_record = nucleotide_index.get(protein_id)
		if nucleotide_record is None:
			missing.append(protein_id)
			continue
		matched.append(nucleotide_record)

	return matched, missing


def make_short_id_map(sequence_ids: list[str]) -> dict[str, str]:
	return {seq_id: f"S{idx:05d}" for idx, seq_id in enumerate(sequence_ids, start=1)}


def _restore_original_ids(value, reverse_id_map: dict[str, str]):
	if isinstance(value, dict):
		restored: dict = {}
		for key, nested_value in value.items():
			new_key = reverse_id_map.get(key, key) if isinstance(key, str) else key
			restored[new_key] = _restore_original_ids(nested_value, reverse_id_map)
		return restored
	if isinstance(value, list):
		return [_restore_original_ids(item, reverse_id_map) for item in value]
	return value


def resolve_yn00_command(command: str) -> str:
	cmd_path = Path(command)
	if cmd_path.is_file():
		return str(cmd_path)

	resolved = shutil.which(command)
	if resolved:
		return resolved

	raise FileNotFoundError(
		"Could not find yn00 executable. Install PAML so 'yn00' is on PATH, "
		"or pass an explicit binary path with --yn00-command /path/to/yn00."
	)


def run_yn00(
	codon_alignment,
	original_ids: list[str],
	yn00_command: str,
	yn00_verbose: bool,
):
	short_id_map = make_short_id_map(original_ids)
	reverse_id_map = {short_id: orig_id for orig_id, short_id in short_id_map.items()}

	for record in codon_alignment:
		record.id = short_id_map[record.id]
		record.name = record.id
		record.description = ""

	tmp_dir = Path(tempfile.mkdtemp(prefix="yn00_"))
	yn00_alignment = tmp_dir / "codon_alignment.phy"
	yn00_out = tmp_dir / "yn00.out"

	try:
		# Write a strict sequential PAML-style file to avoid interleaved PHYLIP parsing issues.
		records = list(codon_alignment)
		if not records:
			raise ValueError("Codon alignment is empty")
		seq_len = len(records[0].seq)
		for rec in records:
			if len(rec.seq) != seq_len:
				raise ValueError("Codon alignment contains sequences with different lengths")

		with yn00_alignment.open("w") as handle:
			handle.write(f" {len(records)} {seq_len}\n")
			for rec in records:
				handle.write(f"{rec.id}\n")
				handle.write(f"{str(rec.seq).upper()}\n")

		runner = yn00.Yn00(
			alignment=str(yn00_alignment),
			working_dir=str(tmp_dir),
			out_file=str(yn00_out),
		)

		parsed_results = runner.run(
			command=yn00_command,
			verbose=yn00_verbose,
			parse=True,
		)
		if parsed_results is None:
			raise ValueError("YN00 returned no parsed results")
		restored_results = _restore_original_ids(dict(parsed_results), reverse_id_map)
	except Exception as exc:
		raw_output = ""
		if yn00_out.exists():
			raw_output = yn00_out.read_text().strip()

		if not yn00_verbose:
			try:
				runner.run(command=yn00_command, verbose=True, parse=False)
			except Exception:
				pass

		details = (
			f"{exc} | Debug files kept in: {tmp_dir}"
		)
		if raw_output:
			details += f" | yn00.out: {raw_output}"
		raise RuntimeError(details) from exc
	else:
		shutil.rmtree(tmp_dir, ignore_errors=True)

	return restored_results


def _fmt(value) -> str:
	if value is None:
		return "NA"
	if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
		return "NA"
	return str(value)


METHODS = ["NG86", "YN00", "LWL85", "LWL85m", "LPB93"]
METHOD_FIELDS: dict[str, list[str]] = {
	"NG86":   ["omega", "dN", "dS"],
	"YN00":   ["S", "N", "t", "kappa", "omega", "dN", "dN SE", "dS", "dS SE"],
	"LWL85":  ["dS", "dN", "w", "S", "N"],
	"LWL85m": ["dS", "dN", "w", "S", "N", "rho"],
	"LPB93":  ["dS", "dN", "w"],
}


def build_tsv_header() -> list[str]:
	cols = ["seq1", "seq2"]
	for method in METHODS:
		for field in METHOD_FIELDS[method]:
			cols.append(f"{method}_{field}")
	return cols


def flatten_results(parsed_results) -> list[dict[str, str]]:
	seen_pairs: set[tuple[str, str]] = set()
	rows: list[dict[str, str]] = []
	for seq1, partners in parsed_results.items():
		if not isinstance(partners, dict):
			continue
		for seq2, methods in partners.items():
			if not isinstance(methods, dict):
				continue
			pair = tuple(sorted((seq1, seq2)))
			if pair in seen_pairs:
				continue
			seen_pairs.add(pair)
			row: dict[str, str] = {"seq1": seq1, "seq2": seq2}
			for method in METHODS:
				method_data = methods.get(method, {})
				for field in METHOD_FIELDS[method]:
					row[f"{method}_{field}"] = _fmt(method_data.get(field))
			rows.append(row)
	return rows


def write_tsv(rows: list[dict[str, str]], output_path: Path) -> None:
	header = build_tsv_header()
	output_path.parent.mkdir(parents=True, exist_ok=True)
	with output_path.open("w", newline="") as handle:
		writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t")
		writer.writeheader()
		writer.writerows(rows)


def main() -> int:
	args = parse_args()

	msa = Path(args.msa)
	fasta = Path(args.fasta)
	output_path = Path(args.output)

	try:
		protein_alignment, protein_ids = read_protein_alignment(msa, args.msa_format)
		nucleotide_index = build_nucleotide_index(fasta)
		matched, missing = find_matching_nucleotides(protein_ids, nucleotide_index)
	except Exception as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		return 1

	if missing:
		print(
			"Missing nucleotide sequence(s) for the following protein ID(s):",
			file=sys.stderr,
		)
		for seq_id in missing:
			print(f"  - {seq_id}", file=sys.stderr)
		return 1

	try:
		codon_alignment = build(protein_alignment, matched)
	except Exception as exc:
		print(f"ERROR: Failed to build codon alignment: {exc}", file=sys.stderr)
		return 1

	try:
		yn00_command = resolve_yn00_command(args.yn00_command)
	except Exception as exc:
		print(f"ERROR: {exc}", file=sys.stderr)
		return 1

	try:
		parsed_results = run_yn00(
			codon_alignment,
			protein_ids,
			yn00_command,
			args.yn00_verbose,
		)
	except Exception as exc:
		print(f"ERROR: Failed to run YN00 pairwise dN/dS: {exc}", file=sys.stderr)
		return 1

	try:
		rows = flatten_results(parsed_results)
		write_tsv(rows, output_path)
	except Exception as exc:
		print(f"ERROR: Failed to write results: {exc}", file=sys.stderr)
		return 1

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
