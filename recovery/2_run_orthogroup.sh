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

# Nextflow puts the pipeline's bin/ on PATH, because ks.py is not in the image.
# Take *only* that one directory. The rest of that line is the host's PATH, and
# exporting it wholesale inside the container shadows the image's own -- which
# is where python3, codeml and yn00 live, under micromamba.
if [[ -z "${KS_BIN:-}" ]]; then
    pathline=$(grep -m1 -oE 'export PATH="[^"]*"' "${SRC}/.command.run" | sed 's/^export PATH="//; s/"$//')
    IFS=':' read -ra entries <<< "${pathline}"
    for e in "${entries[@]}"; do
        [[ -x "${e}/ks.py" ]] && { KS_BIN="${e}"; break; }
    done
fi
[[ -n "${KS_BIN:-}" ]] || { echo "could not find the bin/ directory holding ks.py; set KS_BIN" >&2; exit 1; }

# Reuse the pipeline's container invocation rather than rebuilding the bind
# mounts: the CDS files are symlinks out of the work directory and will not
# resolve inside the container without them.
LAUNCH="${KS_LAUNCH:-$(grep -m1 -oE '(apptainer|singularity) exec [^"]*\.(sif|img)' "${SRC}/.command.run")}"
[[ -n "${LAUNCH}" ]] || { echo "could not derive the container command; set KS_LAUNCH" >&2; exit 1; }
# Nextflow runs this under `env -` so the container starts from the image's own
# environment. Invoked directly it would inherit the host's instead, and the
# host PATH has no python3 -- hence --cleanenv.
[[ "${LAUNCH}" == *--cleanenv* ]] || LAUNCH="${LAUNCH/ exec / exec --cleanenv }"

echo "inputs   : ${SRC}"
echo "bin      : ${KS_BIN}"
echo "cds      : ${KS_CDS}"
echo "args     : ${KS_ARGS} ${KS_ARGS_EXTRA:-}"
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
# Prepend, never replace: \$PATH here is the image's.
export PATH="${KS_BIN}:\$PATH"
# --keep-temp writes under TMPDIR; point it somewhere bound into the container,
# or the kept files vanish with the job.
export TMPDIR="${OUT}"
cd "${OUT}"
# Belt and braces: if the image's PATH did not come through, rebuild it from the
# environment scripts apptainer writes when it converts the Docker image.
if ! command -v python3 >/dev/null 2>&1; then
    set +eu
    for f in /.singularity.d/env/*.sh; do [ -r "\$f" ] && . "\$f"; done
    set -eu
    export PATH="${KS_BIN}:\$PATH"
fi
for prog in python3 codeml yn00 ks.py; do
    command -v "\$prog" >/dev/null || { echo "not on PATH inside the container: \$prog (PATH=\$PATH)" >&2; exit 127; }
done
# \$KS_CDS and \$KS_ARGS word-split deliberately, exactly as the module does.
exec ks.py --msa "${OG}.fas" --tree "${OG}.tree" --paralogs "${OG}.paralogs.txt" \\
    --cds ${KS_CDS} --orthogroup "${OG}" --output "${OG}_ks.tsv" \\
    --threads "${SLURM_CPUS_PER_TASK:-4}" ${KS_ARGS} ${KS_ARGS_EXTRA:-}
INNER

echo "started : $(date)"
time ${LAUNCH} /bin/bash "${OUT}/run.sh"
echo "finished: $(date)"
wc -l "${OUT}/${OG}_ks.tsv"
