#!/usr/bin/env bash
#SBATCH -c 1
#SBATCH --mem=24G
#SBATCH -t 2-00:00:00
#SBATCH -o shadow-%x-%j.out
#SBATCH --constraint=amd9654
#
# Run a second, talkative copy of a tree fit that is already in progress, to
# find out how fast the silent one is going.
#
# The pipeline's control file sets noisy = 0, so a long M0 fit prints nothing
# until it finishes. CODEML under runmode = 0 is deterministic: same binary,
# same alignment, same tree, same control file means the same sequence of
# optimiser rounds, and `noisy` changes only what is printed, not the numerics.
# So a copy started now, on the same CPU, tracks the real job round for round --
# it is simply 59-odd hours behind it.
#
# That gives two things the silent job cannot: minutes per round (so you can
# work out which round the real job is on), and the lnL trace (so you can see
# whether the improvement per round is still decaying geometrically, which
# means convergence is near, or has gone flat, which means it is crawling).
#
# Run it on the SAME constraint as the real job or the timings do not transfer.
#
# Usage:  KS_CODEML=/path/to/codeml sbatch -J OG0000000 5_shadow_fit.sh <run-dir> OG0000000
set -euo pipefail

RUN_DIR="${1:?usage: 5_shadow_fit.sh <run-dir> <OG>}"
OG="${2:?usage: 5_shadow_fit.sh <run-dir> <OG>}"

shopt -s nullglob
src=("${RUN_DIR}"/recovery/"${OG}"/ks_"${OG}"_*/tree)
shopt -u nullglob
[[ ${#src[@]} -gt 0 ]] || {
    echo "No live tree/ directory for ${OG} -- the job has finished or died." >&2
    exit 1
}
SRC="${src[0]}"
OUT="${RUN_DIR}/recovery/${OG}-shadow${KS_SUFFIX:-}"
mkdir -p "${OUT}"

for f in codon_alignment.phy tree.nwk codeml.ctl; do
    [[ -s "${SRC}/${f}" ]] || { echo "missing input: ${SRC}/${f}" >&2; exit 1; }
    cp -f "${SRC}/${f}" "${OUT}/"
done

# noisy = 3 prints a line per optimiser round. Everything else is left exactly
# as the real run has it -- changing any other variable would break the
# round-for-round correspondence this whole approach depends on.
sed -i 's/^\( *noisy *= *\).*/\13/' "${OUT}/codeml.ctl"

# KS_SMALL_DIFF sets minB's e0 -- IF Small_Diff is what Forestry passes there;
# check the minB( call site before relying on this. e0 is tested only as
# `dl < e0 && e <= 0.02`, where dl is the whole round's lnL improvement.
#
# It does NOT change the per-round epsilon. That schedule is internal to minB
# (e /= 2, again if dl < 1, capped at 1e-3 once dl < 0.5) and floors at a
# hardcoded 1e-6, so the expensive rounds happen whatever e0 is set to. What e0
# buys is stopping sooner. Measured on OG0000000: round improvements run 0.130
# (r7), 0.0049 (r8), 0.0020 (r9), 0.0011 (r10), so e0 = 0.01 stops after round
# 8 at ~8 h against ~70-90 h for the full fit, giving up 0.008 lnL out of
# 625249 -- far below what a dS printed to 4 decimals can show.
if [[ -n "${KS_SMALL_DIFF:-}" ]]; then
    sed -i "s/^\\( *Small_Diff *= *\\).*/\\1${KS_SMALL_DIFF}/" "${OUT}/codeml.ctl"
    OUT_NOTE=" (Small_Diff loosened to ${KS_SMALL_DIFF} -- NOT the reference fit)"
fi
grep -E '^ *(noisy|method|runmode|Small_Diff|fix_blength|cleandata|icode) ' "${OUT}/codeml.ctl"
echo "mode     : ${OUT_NOTE:-full convergence}"

CODEML="${KS_CODEML:-codeml}"
command -v "${CODEML}" >/dev/null || { echo "not executable: ${CODEML}" >&2; exit 1; }

echo "node     : $(hostname -s)"
echo "codeml   : $(command -v "${CODEML}")"
echo "inputs   : $(head -1 "${OUT}/codon_alignment.phy")"
echo "started  : $(date '+%F %T')"
echo

cd "${OUT}"
# Unbuffered, timestamped, so each round's wall clock is on the line with it.
stdbuf -oL -eL "${CODEML}" codeml.ctl 2>&1 | stdbuf -oL awk '{ print strftime("%F %T"), $0 }'
