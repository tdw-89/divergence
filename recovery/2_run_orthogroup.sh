#!/usr/bin/env bash
#SBATCH -c 4
#SBATCH --mem=24G
#SBATCH -t 7-00:00:00
#SBATCH -o ks-%x-%j.out
#
# Fit one orthogroup to completion, with no time limit on CODEML.
#
# The M0 tree fit is a single-threaded PAML process whose cost grows as about
# n^1.6 in orthogroup size, so the largest families cannot finish inside a
# shared batch's wall clock however the batches are packed. This runs one on its
# own, reusing the inputs, the CDS list, the arguments and the container image
# the pipeline itself used -- all read back out of the staged work directory, so
# nothing is reconstructed by hand and can drift from what the run did.
#
# Usage:  sbatch -J OG0000000 2_run_orthogroup.sh <run-dir> OG0000000
set -euo pipefail

RUN_DIR="${1:?usage: 2_run_orthogroup.sh <run-dir> <OG>}"
OG="${2:?usage: 2_run_orthogroup.sh <run-dir> <OG>}"
OUT="${RUN_DIR}/recovery/${OG}"

# Any work directory that staged this orthogroup will do; the inputs are
# content-identical wherever they were staged.
SRC=""
while IFS= read -r cand; do
    d=$(dirname "${cand}")
    if [[ -s "${d}/${OG}.paralogs.txt" && -f "${d}/${OG}.tree" && -f "${d}/.command.run" ]]; then
        SRC="${d}"; break
    fi
done < <(find "${RUN_DIR}/work" -maxdepth 3 -name "${OG}.fas")
[[ -n "${SRC}" ]] || { echo "no work directory stages ${OG}" >&2; exit 1; }

# The module exports both of these, so they are the run's own values verbatim.
eval "$(grep -m1 '^export KS_CDS='  "${SRC}/.command.sh")"
eval "$(grep -m1 '^export KS_ARGS=' "${SRC}/.command.sh")"
# Nextflow puts the pipeline's bin/ on PATH; ks.py is not in the image.
BIN_PATH=$(grep -m1 -oE 'export PATH="[^"]*"' "${SRC}/.command.run" || true)
# Reuse the pipeline's container invocation rather than rebuilding the bind
# mounts: the CDS files are symlinks out of the work directory and will not
# resolve inside the container without them.
LAUNCH="${KS_LAUNCH:-$(grep -m1 -oE '(apptainer|singularity) exec [^"]*\.(sif|img)' "${SRC}/.command.run")}"
[[ -n "${LAUNCH}" ]] || { echo "could not derive the container command; set KS_LAUNCH" >&2; exit 1; }

echo "inputs   : ${SRC}"
echo "cds      : ${KS_CDS}"
echo "args     : ${KS_ARGS}"
echo "container: ${LAUNCH}"

mkdir -p "${OUT}"
for f in "${OG}.fas" "${OG}.tree" "${OG}.paralogs.txt" ${KS_CDS}; do
    ln -sfn "$(readlink -f "${SRC}/${f}")" "${OUT}/${f}"
done

# KS_ARGS_EXTRA is appended last so it can override anything the run used --
# e.g. --codeml-command <a locally built codeml>, or
# --keep-temp --skip-pairwise --skip-yn00 to capture just the M0 fit inputs.
cat > "${OUT}/run.sh" <<INNER
set -eu
${BIN_PATH}
cd "${OUT}"
# --keep-temp writes under TMPDIR; point it somewhere bound into the container
# and outside the node's scratch, or the kept files vanish with the job.
export TMPDIR="${OUT}"
# \$KS_CDS and \$KS_ARGS word-split deliberately, exactly as the module does.
exec ks.py --msa "${OG}.fas" --tree "${OG}.tree" --paralogs "${OG}.paralogs.txt" \\
    --cds ${KS_CDS} --orthogroup "${OG}" --output "${OG}_ks.tsv" \\
    --threads "${SLURM_CPUS_PER_TASK:-4}" ${KS_ARGS} ${KS_ARGS_EXTRA:-}
INNER

echo "started : $(date)"
time ${LAUNCH} /bin/bash "${OUT}/run.sh"
echo "finished: $(date)"
wc -l "${OUT}/${OG}_ks.tsv"
