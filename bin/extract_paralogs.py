#!/usr/bin/env python
import argparse
import os
import csv
from Bio import SeqIO

def main():
    parser = argparse.ArgumentParser(description="Extract paralog sequences for a specific species")
    parser.add_argument("--tsv", required=True, help="Path to Orthogroups.tsv")
    parser.add_argument("--seqdir", required=True, help="Path to Orthogroup_Sequences directory")
    parser.add_argument("--species", required=True, help="Target species name (as it appears in OrthoFinder)")
    args = parser.parse_args()

    # Read the Orthogroups.tsv file
    with open(args.tsv, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader)
        
        # Find which column belongs to your target species (prefix match)
        matches = [i for i, col in enumerate(header) if col.startswith(args.species)]
        if not matches:
            raise ValueError(f"Species '{args.species}' not found as a prefix in header: {header}")
        if len(matches) > 1:
            raise ValueError(f"Species '{args.species}' matches multiple columns: {[header[i] for i in matches]}")
        species_idx = matches[0]

        # Process each orthogroup
        for row in reader:
            og_id = row[0]
            # OrthoFinder separates genes with a comma and a space
            species_genes_str = row[species_idx].strip()
            
            if not species_genes_str:
                continue # No genes for this species in this orthogroup
                
            genes = species_genes_str.split(', ')
            
            # We only care about paralogs (2 or more genes in the same species)
            if len(genes) > 1:
                # Open the corresponding sequence file
                seq_file = os.path.join(args.seqdir, f"{og_id}.fa")
                
                if not os.path.exists(seq_file):
                    print(f"Warning: Sequence file {seq_file} not found. Skipping.")
                    continue

                # Parse the sequences and filter for our target genes
                records = list(SeqIO.parse(seq_file, "fasta"))
                paralog_records = [rec for rec in records if rec.id in genes]

                # Write them to a new, species-specific unaligned FASTA
                out_filename = f"{og_id}_{args.species}_unaligned.fasta"
                SeqIO.write(paralog_records, out_filename, "fasta")

if __name__ == "__main__":
    main()