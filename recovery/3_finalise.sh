#!/usr/bin/env bash
# Collect every per-orthogroup result into the run-level codeml/ks.tsv.
#
# This is the only thing the pipeline still owed you: CODEML_BATCH is the last
# process, and the merged table is a `collectFile` -- a concatenation with one
# header, which only materialises at workflow completion. Everything upstream
# already published.
#
# Usage:  3_finalise.sh <run-dir>
set -euo pipefail

RUN_DIR="${1:?usage: 3_finalise.sh <run-dir>}"
DEST="${RUN_DIR}/results/codeml"
mkdir -p "${DEST}/trees"

# Fold in anything produced by 2_run_orthogroup.sh.
shopt -s nullglob
for d in "${RUN_DIR}"/recovery/OG*/; do
    og=$(basename "${d}")
    [[ -s "${d}/${og}_ks.tsv" ]] || { echo "no result yet for ${og}" >&2; continue; }
    cp -f "${d}/${og}_ks.tsv" "${DEST}/"
    for w in dS dN; do
        [[ -s "${d}/${og}_${w}.nwk" ]] && cp -f "${d}/${og}_${w}.nwk" "${DEST}/trees/"
    done
done

parts=("${DEST}"/*_ks.tsv)
(( ${#parts[@]} )) || { echo "no *_ks.tsv under ${DEST}" >&2; exit 1; }

# One header, then every data row. Guard against a part whose schema drifted.
header=$(head -1 "${parts[0]}")
tmp=$(mktemp)
printf '%s\n' "${header}" > "${tmp}"
bad=0
for f in "${parts[@]}"; do
    if [[ "$(head -1 "${f}")" != "${header}" ]]; then
        echo "SKIPPED (header differs): ${f}" >&2; bad=$((bad + 1)); continue
    fi
    tail -n +2 "${f}" >> "${tmp}"
done
mv "${tmp}" "${DEST}/ks.tsv"

echo "merged ${#parts[@]} orthogroup table(s), ${bad} skipped"
echo "rows: $(( $(wc -l < "${DEST}/ks.tsv") - 1 ))"
echo "orthogroups represented: $(tail -n +2 "${DEST}/ks.tsv" | cut -f1 | sort -u | wc -l)"
