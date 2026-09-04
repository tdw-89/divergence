#!/usr/bin/env python3
"""Profile an orthogroup alignment before deciding to carve it up.

Reads the protein alignment MAFFT_BATCH produced (<OG>.fas) and, optionally,
the manifest extract_paralogs.py wrote beside it. Pure stdlib: this is meant to
run anywhere, including a login node with no biopython.

The question it answers is which of two very different things is making a
CODEML fit expensive:

  * sequence count -- the M0 fit is ~O(n^1.6) in sequences, so splitting a
    family pays badly (a 2-way split saves 34%, and you need a 16-way split to
    save 81%), and every cut severs the reported pairs that span it;
  * alignment width -- cost is LINEAR in columns, so trimming columns that
    almost nothing occupies is a far better trade, loses no sequences, and
    breaks no pairs.

Pair coverage is deliberately measured as shared / min(len_a, len_b), never
over the alignment width. Dividing by the width makes a healthy family of
variable-length proteins look identical to a broken one; that error once
produced a false "MCL-chained orthogroups" diagnosis in this project.
"""
import argparse, random, statistics as st, sys
from collections import Counter


def read_fasta(path):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            elif name is not None:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf)
    if not seqs:
        sys.exit(f"no records in {path}")
    return seqs


def quantiles(xs):
    xs = sorted(xs)
    q = lambda f: xs[min(len(xs) - 1, int(f * len(xs)))]
    return xs[0], q(0.25), q(0.50), q(0.75), xs[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("alignment")
    ap.add_argument("--manifest", help="<OG>.paralogs.txt, to count reportable pairs")
    ap.add_argument("--pairs", type=int, default=4000, help="pairs to sample for coverage")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    seqs = read_fasta(a.alignment)
    ids = list(seqs)
    width = len(seqs[ids[0]])
    bad = [i for i in ids if len(seqs[i]) != width]
    if bad:
        sys.exit(f"not an alignment: {len(bad)} records differ in length, e.g. {bad[0]}")

    masks = {i: [c != "-" for c in s] for i, s in seqs.items()}
    lens = {i: sum(m) for i, m in masks.items()}
    total = sum(lens.values())

    print(f"file          : {a.alignment}")
    print(f"sequences     : {len(ids)}")
    print(f"columns       : {width}")
    print(f"gap fraction  : {1 - total / (len(ids) * width):.4f}")
    lo, q1, med, q3, hi = quantiles(lens.values())
    print(f"ungapped len  : min {lo}  q1 {q1}  median {med}  q3 {q3}  max {hi}")
    print(f"width / median: {width / med:.1f}x"
          f"   <- a healthy family sits near 1-2x")

    occ = [0] * width
    for m in masks.values():
        for c, present in enumerate(m):
            if present:
                occ[c] += 1

    print("\ncolumn occupancy -- what trimming would buy")
    print("  min seqs   columns kept   width    predicted fit cost (linear in width)")
    for thr_pct in (0, 1, 2, 5, 10, 20, 50):
        thr = max(1, len(ids) * thr_pct // 100)
        kept = sum(1 for o in occ if o >= thr)
        print(f"  >={thr_pct:>3}% ({thr:>4})   {kept:>8}   {kept / width:>6.1%}"
              f"   {kept / width:>6.2f}x")

    rng = random.Random(a.seed)
    sample = ([(x, y) for n, x in enumerate(ids) for y in ids[n + 1:]]
              if len(ids) * (len(ids) - 1) // 2 <= a.pairs
              else [tuple(rng.sample(ids, 2)) for _ in range(a.pairs)])
    covs = []
    for x, y in sample:
        mx, my = masks[x], masks[y]
        shared = sum(1 for c in range(width) if mx[c] and my[c])
        denom = min(lens[x], lens[y])
        if denom:
            covs.append(shared / denom)
    if covs:
        lo, q1, med, q3, hi = quantiles(covs)
        print(f"\npair coverage (shared / min length), {len(covs)} pairs sampled")
        print(f"  min {lo:.3f}  q1 {q1:.3f}  median {med:.3f}  q3 {q3:.3f}  max {hi:.3f}")
        for t in (0.5, 0.2, 0.05):
            n = sum(1 for c in covs if c < t)
            print(f"  below {t:<5}: {n:>6} ({n / len(covs):.1%})")

    if a.manifest:
        by_species = Counter()
        with open(a.manifest) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    by_species[parts[0]] += 1
        print("\nreportable pairs by target species (manifest)")
        tot = 0
        for sp, n in sorted(by_species.items(), key=lambda kv: -kv[1]):
            p = n * (n - 1) // 2
            tot += p
            print(f"  {sp:<40} {n:>5} genes  {p:>8} pairs")
        print(f"  {'TOTAL':<40} {'':>5}         {tot:>8} pairs")


if __name__ == "__main__":
    main()
