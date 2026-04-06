#!/usr/bin/env python3

import argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

codon_map = {
    'A': ['GCT', 'GCC', 'GCA', 'GCG'],
    'C': ['TGT', 'TGC'],
    'D': ['GAT', 'GAC'],
    'E': ['GAA', 'GAG'],
    'F': ['TTT', 'TTC'],
    'G': ['GGT', 'GGC', 'GGA', 'GGG'],
    'H': ['CAT', 'CAC'],
    'I': ['ATT', 'ATC', 'ATA'],
    'K': ['AAA', 'AAG'],
    'L': ['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'],
    'M': ['ATG'],
    'N': ['AAT', 'AAC'],
    'P': ['CCT', 'CCC', 'CCA', 'CCG'],
    'Q': ['CAA', 'CAG'],
    'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],
    'T': ['ACT', 'ACC', 'ACA', 'ACG'],
    'V': ['GTT', 'GTC', 'GTA', 'GTG'],
    'W': ['TGG'],
    'Y': ['TAT', 'TAC'],
    '*': ['TAA', 'TAG', 'TGA'],
    'X': ['NNN'] 
}

def generate_plausible_cds(input_fasta, output_fasta):
    records = []
    for record in SeqIO.parse(input_fasta, "fasta"):
        clean_seq = str(record.seq).strip()
        
        nt_seq = "".join([random.choice(codon_map.get(aa.upper(), ['NNN'])) for aa in clean_seq])
        
        new_record = SeqRecord(Seq(nt_seq), id=record.id, description=record.description)
        records.append(new_record)
        
    SeqIO.write(records, output_fasta, "fasta")
    print(f"Successfully generated plausible CDS FASTA at {output_fasta}")

def main():
    parser = argparse.ArgumentParser(description="Generate a plausible nucleotide CDS FASTA from an unaligned protein FASTA.")
    parser.add_argument("-i", "--input", required=True, help="Input unaligned protein FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output nucleotide CDS FASTA file")
    args = parser.parse_args()

    generate_plausible_cds(args.input, args.output)

if __name__ == "__main__":
    main()