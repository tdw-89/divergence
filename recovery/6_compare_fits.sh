#!/usr/bin/env bash
# Compare two ks.py runs of the same orthogroup that differ only in codeml's
# convergence tolerance, and report whether the difference is visible in what
# the pipeline actually publishes.
#
# lnL is the wrong yardstick for this. CLAUDE.md already records a case where a
# 0.11 lnL gap moved 10 of 117 branches' dS and collapsed one to zero, so a
# small lnL delta is not on its own evidence that the estimates agree. What
# matters is tree_dS per reported pair, which is a path sum over the fitted
# branch table and is printed by codeml to 4 decimals.
#
#   recovery/6_compare_fits.sh <a_ks.tsv> <b_ks.tsv> [tolerance]
#
# Default tolerance 5e-5 is half of codeml's own last printed digit: differences
# below it cannot be represented in the branch table the sums are built from.
set -uo pipefail

A="${1:?usage: 6_compare_fits.sh <a_ks.tsv> <b_ks.tsv> [tolerance]}"
B="${2:?usage: 6_compare_fits.sh <a_ks.tsv> <b_ks.tsv> [tolerance]}"
TOL="${3:-5e-5}"

for f in "${A}" "${B}"; do
    [[ -s "${f}" ]] || { echo "missing or empty: ${f}" >&2; exit 1; }
done

python3 - "${A}" "${B}" "${TOL}" <<'PY'
import csv, sys

a_path, b_path, tol = sys.argv[1], sys.argv[2], float(sys.argv[3])
KEY = ("orthogroup", "species", "gene_a", "gene_b")
# m0_* are per-orthogroup, not per-pair; report them separately.
COLS = ("tree_dS", "tree_dN", "tree_omega")

def load(path):
    with open(path, newline="") as fh:
        rows = {tuple(r[k] for k in KEY): r for r in csv.DictReader(fh, delimiter="\t")}
    if not rows:
        sys.exit(f"no data rows in {path}")
    return rows

a, b = load(a_path), load(b_path)
only_a, only_b = sorted(a.keys() - b.keys()), sorted(b.keys() - a.keys())
shared = sorted(a.keys() & b.keys())

print(f"pairs: {len(a)} in A, {len(b)} in B, {len(shared)} shared")
for label, missing in (("A only", only_a), ("B only", only_b)):
    if missing:
        print(f"  !! {len(missing)} {label} -- the two runs do not describe the same pair set")
        for k in missing[:5]:
            print(f"     {'/'.join(k)}")

for col in COLS:
    diffs, na_mismatch = [], 0
    for k in shared:
        va, vb = a[k][col], b[k][col]
        if va == "NA" or vb == "NA":
            na_mismatch += va != vb
            continue
        diffs.append((abs(float(va) - float(vb)), k, float(va), float(vb)))
    print(f"\n{col}: {len(diffs)} comparable, {na_mismatch} NA-mismatched")
    if not diffs:
        continue
    diffs.sort(reverse=True)
    over = [d for d in diffs if d[0] > tol]
    vals = [d[0] for d in diffs]
    print(f"  max {vals[0]:.3e}   median {sorted(vals)[len(vals)//2]:.3e}"
          f"   over {tol:g}: {len(over)}/{len(diffs)}")
    for d, k, va, vb in diffs[:3]:
        rel = f"{d / abs(va) * 100:.2f}%" if va else "n/a"
        print(f"    {'/'.join(k[1:])}: {va:.6f} vs {vb:.6f}  (d={d:.3e}, {rel})")

for col in ("m0_lnL", "m0_kappa", "m0_omega"):
    va, vb = a[shared[0]][col], b[shared[0]][col]
    print(f"{col:<10} A={va:<20} B={vb}")
PY
