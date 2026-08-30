#!/usr/bin/env bash
# Recover the completed orthogroups from CODEML_BATCH tasks that SLURM killed.
#
# A task that exits non-zero never runs publishDir, so its per-orthogroup
# results stay in the work directory and never reach results/. But CODEML_BATCH
# writes <OG>_ks.tsv as each orthogroup finishes, so a batch killed part way
# through still holds complete, valid results for everything that did finish --
# in this run, 19 of every 20.
#
# Usage:  1_harvest.sh <run-dir>          # e.g. .../divergence/fish
set -euo pipefail

RUN_DIR="${1:?usage: 1_harvest.sh <run-dir>}"
WORK="${RUN_DIR}/work"
DEST="${RUN_DIR}/results/codeml"
mkdir -p "${DEST}/trees"

harvested=0
skipped=0
missing=()

# Every CODEML_BATCH work directory that exited non-zero.
while IFS= read -r exitfile; do
    d=$(dirname "${exitfile}")
    [[ -f "${d}/.command.sh" ]] || continue
    grep -q 'ks.py' "${d}/.command.sh" || continue
    [[ "$(cat "${exitfile}")" == "0" ]] && continue

    for aln in "${d}"/*.fas; do
        [[ -e "${aln}" ]] || continue
        og=$(basename "${aln}"); og="${og%%.*}"
        tsv="${d}/${og}_ks.tsv"

        if [[ ! -s "${tsv}" ]]; then
            missing+=("${og}")
            continue
        fi
        if [[ -e "${DEST}/${og}_ks.tsv" ]]; then
            skipped=$((skipped + 1))
            continue
        fi
        # A file written by a completed ks.py call has a header and at least
        # one data row, and every row has the header's column count.
        ncol=$(head -1 "${tsv}" | awk -F'\t' '{print NF}')
        if ! awk -F'\t' -v n="${ncol}" 'NF != n {bad=1} END {exit (bad || NR < 2)}' "${tsv}"; then
            echo "SUSPECT (not harvested): ${tsv}" >&2
            continue
        fi

        cp -n "${tsv}" "${DEST}/"
        for w in dS dN; do
            [[ -s "${d}/${og}_${w}.nwk" ]] && cp -n "${d}/${og}_${w}.nwk" "${DEST}/trees/"
        done
        harvested=$((harvested + 1))
    done
done < <(find "${WORK}" -maxdepth 3 -name .exitcode)

echo "harvested ${harvested} orthogroup(s); ${skipped} already published"
if ((${#missing[@]})); then
    printf 'still missing (%d):\n' "${#missing[@]}"
    printf '  %s\n' "${missing[@]}" | sort -u
fi
