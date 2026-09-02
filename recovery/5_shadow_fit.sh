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
OUT="${RUN_DIR}/recovery/${OG}-shadow"
mkdir -p "${OUT}"

for f in codon_alignment.phy tree.nwk codeml.ctl; do
    [[ -s "${SRC}/${f}" ]] || { echo "missing input: ${SRC}/${f}" >&2; exit 1; }
    cp -f "${SRC}/${f}" "${OUT}/"
done

# noisy = 3 prints a line per optimiser round. Everything else is left exactly
# as the real run has it -- changing any other variable would break the
# round-for-round correspondence this whole approach depends on.
sed -i 's/^\( *noisy *= *\).*/\13/' "${OUT}/codeml.ctl"
grep -E '^ *(noisy|method|runmode|Small_Diff|fix_blength|cleandata|icode) ' "${OUT}/codeml.ctl"

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
