#!/usr/bin/env python3
"""Simulate a small gene family with known synonymous divergence.

Mutations are applied only at synonymous positions, so the true omega is ~0 and
true dS along a path is roughly (mutations applied) / (synonymous sites).
Topology:  ((A,B):inner, (C,D):inner, OUT)  with A,B,C,D the 'target paralogs'.
"""
import random
from pathlib import Path
from Bio.Data import CodonTable
from Bio.Seq import Seq

random.seed(20260821)
TABLE = CodonTable.unambiguous_dna_by_id[1]

# codon -> list of synonymous alternatives
SYN = {}
for codon, aa in TABLE.forward_table.items():
    SYN.setdefault(aa, []).append(codon)
SYNONYMS = {c: [x for x in SYN[aa] if x != c] for c, aa in TABLE.forward_table.items()}

NCOD = 400
ancestor = [random.choice([c for c in TABLE.forward_table if SYNONYMS[c]]) for _ in range(NCOD)]

def mutate(seq, n):
    seq = list(seq)
    positions = random.sample(range(len(seq)), n)
    for pos in positions:
        alts = SYNONYMS.get(seq[pos])
        if alts:
            seq[pos] = random.choice(alts)
    return seq

# branch lengths expressed as number of synonymous codon swaps
root      = ancestor
node_ab   = mutate(root, 20)
node_cd   = mutate(root, 20)
seqs = {
    "geneA": mutate(node_ab, 10),
    "geneB": mutate(node_ab, 10),
    "geneC": mutate(node_cd, 10),
    "geneD": mutate(node_cd, 10),
    "outgroup1": mutate(root, 40),
    "outgroup2": mutate(root, 45),
}

out = Path(__file__).parent
with (out / "cds_target.fa").open("w") as fh, (out / "cds_other.fa").open("w") as oh:
    for name, codons in seqs.items():
        handle = oh if name.startswith("outgroup") else fh
        handle.write(f">{name}\n{''.join(codons)}\n")

# protein 'alignment': no indels were simulated, so translations line up already
with (out / "aln.faa").open("w") as fh:
    for name, codons in seqs.items():
        fh.write(f">{name}\n{Seq(''.join(codons)).translate()}\n")

(out / "tree.nwk").write_text(
    "(((geneA,geneB),(geneC,geneD)),outgroup1,outgroup2);\n"
)
(out / "paralogs.txt").write_text("geneA\ngeneB\ngeneC\ngeneD\n")

nsyn = NCOD * 1.0  # rough: ~1 synonymous site per codon at 3rd position
print(f"{NCOD} codons; expected dS A-B ~ {20/ (NCOD*0.85):.4f}, A-C ~ {60/(NCOD*0.85):.4f}")
