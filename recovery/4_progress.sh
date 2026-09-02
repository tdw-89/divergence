#!/usr/bin/env bash
# Ask a still-running 2_run_orthogroup.sh job how far through it is.
#
# ks.py sets TMPDIR to the job's own output directory, so CODEML's live working
# directory sits at <run-dir>/recovery/<OG>/ks_<OG>_*/tree/ for as long as the
# fit is running -- it is deleted on exit unless --keep-temp was given. That
# directory is the only progress signal there is: CODEML's console output is
# captured by a pipe in _run_paml and the control file sets noisy = 0, so
# nothing reaches the SLURM log until the whole orthogroup finishes.
#
# The file that matters is `rub`, which PAML appends to as the optimiser runs.
# Under method = 1 each round costs about the same -- one line search per branch
# -- so rounds/hour is near-constant and extrapolating over it is legitimate,
# which is not true of the method = 0 optimiser.
#
# Usage:  4_progress.sh <run-dir> <OG> [sample-seconds]
set -euo pipefail

RUN_DIR="${1:?usage: 4_progress.sh <run-dir> <OG> [sample-seconds]}"
OG="${2:?usage: 4_progress.sh <run-dir> <OG> [sample-seconds]}"
INTERVAL="${3:-300}"

OUT="${RUN_DIR}/recovery/${OG}"
[[ -d "${OUT}" ]] || { echo "no such job directory: ${OUT}" >&2; exit 1; }

shopt -s nullglob
roots=("${OUT}"/ks_"${OG}"_*)
shopt -u nullglob
if [[ ${#roots[@]} -eq 0 ]]; then
    echo "No live ks_${OG}_* directory under ${OUT}."
    echo "Either the job has finished (check for ${OG}_ks.tsv) or it died before"
    echo "the codon alignment was built."
    ls -l "${OUT}" || true
    exit 0
fi
ROOT="${roots[0]}"
TREE="${ROOT}/tree"

echo "=== job ==="
echo "work root : ${ROOT}"
echo "started   : $(date -r "${ROOT}" '+%F %T')  ($(( ( $(date +%s) - $(stat -c %Y "${ROOT}") ) / 3600 ))h ago)"
squeue -h -o '%i %T %M %L %N' -n "${OG}" 2>/dev/null || true

echo
echo "=== phase ==="
shopt -s nullglob
pairdirs=("${ROOT}"/pair_*)
shopt -u nullglob
if [[ ${#pairdirs[@]} -gt 0 ]]; then
    echo "PAST the tree fit -- ${#pairdirs[@]} pairwise/YN00 directories exist."
    echo "Those run at ~0.088 s/pair across --threads workers, so the rest is minutes."
    exit 0
fi
if [[ ! -d "${TREE}" ]]; then
    echo "tree/ does not exist yet: still building the codon alignment."
    exit 0
fi
echo "IN the M0 tree fit (no pair_* directories yet)."

echo
echo "=== problem size ==="
if [[ -s "${TREE}/codon_alignment.phy" ]]; then
    read -r nseq nsite < <(head -1 "${TREE}/codon_alignment.phy")
    echo "sequences : ${nseq}"
    echo "codons    : $(( nsite / 3 ))"
fi
[[ -s "${TREE}/lnf" ]] && echo "lnf lines : $(wc -l < "${TREE}/lnf")  (roughly the site-pattern count)"

echo
echo "=== files ==="
ls -lt --time-style='+%F %T' "${TREE}"

echo
echo "=== rub (optimiser trace) ==="
if [[ -s "${TREE}/rub" ]]; then
    echo "-- first 3 --"; head -3 "${TREE}/rub"
    echo "-- last 5 --";  tail -5 "${TREE}/rub"
    n1=$(wc -l < "${TREE}/rub")
    echo
    echo "sampling for ${INTERVAL}s ..."
    sleep "${INTERVAL}"
    n2=$(wc -l < "${TREE}/rub")
    echo "lines: ${n1} -> ${n2}  (+$(( n2 - n1 )) in ${INTERVAL}s)"
    if [[ $(( n2 - n1 )) -gt 0 ]]; then
        elapsed=$(( $(date +%s) - $(stat -c %Y "${ROOT}") ))
        awk -v a="${n1}" -v b="${n2}" -v s="${INTERVAL}" -v e="${elapsed}" -v n="${n2}" \
            'BEGIN {
                 rate = (b - a) / s * 3600;
                 printf "rate     : %.2f rounds/hour (%.1f min per round)\n", rate, 60 / rate;
                 printf "elapsed  : %.1f h for %d rounds (%.1f min/round averaged)\n", e/3600, n, e/60/n;
                 printf "budget   : %.1f h left of the 7-day limit = room for ~%.0f more rounds\n", \
                        (604800 - e)/3600, (604800 - e)/3600 * rate;
             }'
    else
        echo "NO PROGRESS in ${INTERVAL}s -- either rub is not the live trace in this"
        echo "PAML build, or a single round is longer than the sample window."
    fi
else
    echo "rub is absent or empty."
    echo "Check whether this build writes it:  grep -n 'frub' /path/to/paml/src/*.c"
fi

echo
echo "=== mlc growth ==="
if [[ -e "${TREE}/mlc" ]]; then
    stat -c 'size %s bytes, modified %y' "${TREE}/mlc"
else
    echo "mlc not created yet"
fi
