# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`divergence` is a Nextflow (DSL2) pipeline built from the nf-core template (`nf-core/tools` 3.5.2, `is_nfcore: false`, template org `uml`). It takes protein FASTAs for several species, runs OrthoFinder to build orthogroups, pulls out the *paralogs of one target species*, aligns each paralog set, and estimates pairwise dN/dS with PAML's YN00.

## Commands

Run the pipeline (test profile uses four Mycoplasma proteomes in `assets/test_fastas/`):

```bash
nextflow run . -profile test,docker --outdir results
```

Real run — all three params are required, plus `--outdir`:

```bash
nextflow run . -profile docker \
  --input samplesheet.csv \
  --target_species <species_name> \
  --nuc_ref <target_species_cds.fa> \
  --outdir results
```

Pipeline-level tests (nf-test, `testsDir "."`, profile `test`, nf-core module/subworkflow tests ignored):

```bash
nf-test test                          # everything
nf-test test tests/default.nf.test    # the single pipeline test
nf-test test --profile test,docker    # override container engine
nf-test test --update-snapshot        # regenerate tests/default.nf.test.snap
```

Python unit tests for `bin/dnds.py` (stdlib `unittest`, needs biopython):

```bash
python -m unittest discover -s bin/tests   # or: pytest bin/tests
```

Linting: `pre-commit run --all-files` (prettier + whitespace fixers) and `nf-core pipelines lint`. CI runs nf-test across `{conda, docker, singularity} × {25.04.0, latest-everything}`, sharded.

## Architecture

`main.nf` → `PIPELINE_INITIALISATION` (nf-schema validation of `assets/schema_input.json`, producing `[meta, fasta]`) → `workflows/divergence.nf` (the `DIVERGENCE` workflow, where all real logic lives) → `PIPELINE_COMPLETION`.

Inside `DIVERGENCE` the flow is strictly linear (see `diag.md` for the mermaid version):

1. **ORTHOFINDER** (`modules/nf-core/orthofinder`) — every samplesheet FASTA is `collect()`ed into one channel item tagged `[id: 'orthofinder']`; the prior-run input is always empty.
2. **EXTRACT_PARALOGS** (local) — runs `bin/extract_paralogs.py` against `Orthogroups/Orthogroups.tsv` + `Orthogroup_Sequences/` from the OrthoFinder output dir. It picks the header column whose name *starts with* `--target_species` (errors on zero or multiple matches), keeps only orthogroups where that species has ≥2 genes, and emits one `<OG>_<species>_unaligned.fasta` per such orthogroup. Output has no meta map — it is a bare `path` channel.
3. **MAFFT_BATCH** (local) — `.flatten().collate(params.batch_size)` groups the per-orthogroup FASTAs into batches; the process loops over the batch in bash, one `mafft` call per file, writing `<prefix>.fas`. `--anysymbol` is set via `ext.args` in `conf/modules.config`.
4. **DNDS_BATCH** (local) — same flatten/collate batching over the alignments; loops calling `bin/dnds.py` once per MSA with the shared `--nuc_ref` nucleotide FASTA, emitting `<prefix>_dnds.tsv`.

**Batching is the key structural idea**: orthogroup counts run into the thousands, so `params.batch_size` (default 20) controls how many files each MAFFT/dN/dS task handles inside a bash loop, rather than spawning one Nextflow task per orthogroup. This means the modules take plain `path(files)` collections instead of `tuple val(meta), path(file)`, and per-file failures must be handled inside the shell loop.

**Failure semantics in DNDS_BATCH** are deliberate and were arrived at over several commits: a single MSA that fails or produces an empty TSV logs a warning and continues; the task only exits non-zero if *every* MSA in the batch failed. `tsv` is `optional: true`. Don't tighten this to fail on individual MSAs.

`bin/dnds.py` (the only non-trivial script):
- Maps protein MSA IDs to nucleotide records by **exact `record.id` match** against `--nuc_ref`. Strict matching is intentional — a previous attempt at fuzzy ID parsing was reverted. IDs with no nucleotide match are dropped with a warning; fewer than 2 remaining sequences is an error.
- Builds a codon alignment with `Bio.codonalign.build`, renames sequences to short `S00001`-style IDs (PAML name-length limits) and restores the originals afterwards.
- Writes a strict sequential PAML-style `.phy` by hand rather than using Biopython's PHYLIP writer.
- If Biopython's YN00 parser raises `IndexError` on the whole-alignment run (usually too few informative sites), it falls back to `_run_yn00_pairwise`, running YN00 on each pair in its own temp dir so only the affected pairs come out as `NA`.
- Adds bidirectional pairwise % identity columns computed from the codon alignment alongside the YN00 fields.
- Note: this file is **tab-indented**; every other script in `bin/` uses 4 spaces. Match whatever file you are editing.

`bin/generate_cds.py`, `bin/longest_isoform.py`, and `bin/get_matching_cds.py` are standalone data-prep helpers, not called by any process — they exist to build the `--nuc_ref` CDS FASTA and to filter proteomes to one isoform per gene before running the pipeline. `get_matching_cds.py` accepts a CDS record only if the protein ID appears in the defline *and* the CDS translates exactly to the protein sequence.

## Containers

`MAFFT_BATCH` and `EXTRACT_PARALOGS` use standard biocontainers. `DNDS_BATCH` uses a custom image, `docker.io/tdw0student0uml/dnds:1.0`, built from `modules/local/dnds_batch/Dockerfile` (micromamba + biopython 1.84 + paml 4.10.7). Changing the dN/dS dependencies means rebuilding and pushing that image, not just editing the conda directive.

## Caveats

- `README.md` and `docs/usage.md` are still unmodified nf-core template boilerplate — the fastq/samplesheet examples there describe a read-based pipeline and do **not** reflect this one. The real input format is `species,fasta` (see `assets/test_samplesheet.csv` and `assets/schema_input.json`). Don't treat those docs as a source of truth; updating them is outstanding work.
- `.nf-core.yml` disables a large set of template lint checks (multiqc, igenomes, license, changelog, logos). Don't add the corresponding files back to satisfy a lint error.
- `null/`, `work/`, `diagnose/`, and `tests/out/` are gitignored run detritus.
