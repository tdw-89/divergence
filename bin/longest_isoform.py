#!/usr/bin/env python3
import argparse
import os
import re
import sys
import warnings
from Bio import SeqIO

def main():
    # Set up argument parsing
    parser = argparse.ArgumentParser(
        description="Filter a protein FASTA file to keep only the longest sequence per gene."
    )
    parser.add_argument(
        "input_fasta",
        help="Path to the mandatory input FASTA file."
    )
    parser.add_argument(
        "-o", "--output",
        help="Path to the output FASTA file. (Defaults to <input_path>_longest.<ext>)"
    )
    parser.add_argument(
        "-r", "--regex",
        default=r"gene:[A-Z0-9\.]+",
        help="Regular expression pattern to identify the gene ID in the defline. Default: 'gene:[A-Z0-9\\.]+'"
    )

    args = parser.parse_args()

    # 1. Validate the regular expression
    try:
        gene_pattern = re.compile(args.regex)
    except re.error as e:
        raise ValueError(f"Invalid regular expression provided: '{args.regex}'. Error: {e}")

    # 2. Determine the output path if not explicitly provided
    if args.output:
        output_fasta = args.output
    else:
        base, ext = os.path.splitext(args.input_fasta)
        output_fasta = f"{base}_longest{ext}"

    longest_records = {}

    # 3. Load the protein FASTA and iterate through records
    print(f"Parsing '{args.input_fasta}'...")
    try:
        for record in SeqIO.parse(args.input_fasta, "fasta"):
            # Search for the gene ID in the defline (record.description)
            match = gene_pattern.search(record.description)
            
            # If the pattern doesn't match, warn the user and skip the record
            if not match:
                warnings.warn(
                    f"\nRegex pattern '{args.regex}' did not match anything in defline: '{record.description}'. Skipping record."
                )
                continue
            
            gene_id = match.group(0)
            
            # 4. Add to dictionary or replace if the current sequence is longer
            if gene_id not in longest_records:
                longest_records[gene_id] = record
            else:
                if len(record.seq) > len(longest_records[gene_id].seq):
                    longest_records[gene_id] = record
                    
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_fasta}' not found.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"An error occurred while parsing the FASTA file: {e}", file=sys.stderr)
        sys.exit(1)

    # 5. Write the filtered records to the new FASTA file
    try:
        SeqIO.write(longest_records.values(), output_fasta, "fasta")
        print(f"Success! Wrote {len(longest_records)} unique longest protein sequences to '{output_fasta}'.")
    except Exception as e:
        print(f"An error occurred while writing to the output file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
