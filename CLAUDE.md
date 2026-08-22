# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

`divergence` is a Nextflow (DSL2) pipeline built from the nf-core template (`nf-core/tools` 3.5.2, `is_nfcore: false`, template org `uml`). It takes protein FASTAs and matching CDS FASTAs for several species, runs OrthoFinder, selects the orthogroups containing paralogs of one or more target species, aligns them, and estimates **pairwise dS (Ks) between those paralogs** with PAML's CODEML.

The scientific goal is _time since duplication_, with dS as the proxy. The pipeline deliberately stops there: it does not compute last-common-ancestor brackets or ohnolog calls, but it does publish the OrthoFinder outputs those analyses need (`Gene_Duplication_Events/`, `Phylogenetic_Hierarchical_Orthogroups/`, `Resolved_Gene_Trees/`). Synteny and LCA analysis are done downstream, outside the pipeline.

## Commands

Run the pipeline (test profile uses four Mycoplasma proteomes in `assets/test_fastas/`):

```bash
nextflow run . -profile test,docker --outdir results
```

Real run — both params are required, plus `--outdir`:

```bash
nextflow run . -profile docker \
  --input samplesheet.csv \
  --target_species <species_name> \
  --outdir results
```

Samplesheet is `species,fasta,cds` (see `assets/test_samplesheet.csv`). **Every species needs a CDS file**, not just the target: codon alignments are built for whole orthogroups. CDS record IDs must match the protein IDs exactly — `bin/get_matching_cds.py` produces files in that form.

Pipeline-level tests (nf-test, `testsDir "."`, profile `test`, nf-core module/subworkflow tests ignored):

```bash
nf-test test                                    # everything
nf-test test tests/default.nf.test              # the pipeline test
nf-test test tests/modules/local/codeml_batch   # the CODEML_BATCH module test
nf-test test --profile test,docker              # override container engine
```

Neither test uses a snapshot: the numbers are ML estimates and are not bit-stable
across PAML builds or platforms, so they assert column schema, counts and value
ranges instead. There is no `.snap` file to regenerate.

**nf-test needs Java 17–21** — it fails with `Unsupported class file major version`
on newer JDKs, even ones Nextflow itself accepts. Nextflow needs 17–24. If the
default JDK is newer than either, set `JAVA_HOME` accordingly.

Python unit tests (stdlib `unittest`, needs biopython):

```bash
python -m unittest discover -s bin/tests   # or: pytest bin/tests
```

There is no local biopython/PAML environment on the dev machine; the easiest way to run the Python tests, or to exercise `ks.py` by hand, is inside the pipeline's own image:

```bash
docker run --rm --platform linux/amd64 -v "$PWD:/w" -w /w \
  docker.io/tdw0student0uml/dnds:1.0 python3 -m unittest discover -s bin/tests
```

Linting: `pre-commit run --all-files` (prettier + whitespace fixers) and `nf-core pipelines lint`. CI runs nf-test across `{conda, docker, singularity} × {25.04.0, latest-everything}`, sharded.

## Architecture

`main.nf` → `PIPELINE_INITIALISATION` (nf-schema validation of `assets/schema_input.json`, producing `[meta, fasta, cds]`) → `workflows/divergence.nf` (the `DIVERGENCE` workflow, where all real logic lives) → `PIPELINE_COMPLETION`.

1. **ORTHOFINDER** (`modules/nf-core/orthofinder`) — every samplesheet proteome is `collect()`ed into one channel item tagged `[id: 'orthofinder']`; the prior-run input is always empty. The whole output directory is published, which is what makes downstream LCA/ohnolog analysis possible.
2. **EXTRACT_PARALOGS** (local) — runs `bin/extract_paralogs.py`. `--target_species` is a comma-separated list or `all`; each name is matched as a header-column _prefix_ (errors on zero or multiple matches). An orthogroup is kept when any target species has ≥ `--min_paralogs` genes there, and `--max_target_paralogs` drops an over-large species from an orthogroup rather than discarding the orthogroup. For each, emits three files sharing an `<OG>` stem: `<OG>.faa` (the **whole** orthogroup, all species), `<OG>.tree` (gene tree), `<OG>.paralogs.txt` (a `species<TAB>gene` manifest). Also writes `orthogroup_summary.tsv` covering every orthogroup considered, kept or not.
3. **MAFFT_BATCH** (local) — `.flatten().collate(params.batch_size)` batches the FASTAs; the process loops in bash, one `mafft` call per file, writing `<OG>.fas`. `--anysymbol` via `ext.args`.
4. **CODEML_BATCH** (local) — loops calling `bin/ks.py` once per orthogroup, emitting `<OG>_ks.tsv` plus `<OG>_dS.nwk` / `<OG>_dN.nwk` (published under `codeml/trees/`).

**Why whole orthogroups.** The extra species are never reported on. They exist to break up long branches: a paralog pair whose duplication predates several speciations sits at the end of one enormous branch, where dS is saturated and the multiple-hit correction has nothing to work with. Interposing orthologs subdivides that path into individually short branches. Do not "optimise" this back to target-species-only alignments — it is the point of the design.

**Batching** — orthogroup counts run into the thousands, so `params.batch_size` (default 20) controls how many orthogroups each MAFFT/CODEML task handles in a bash loop, rather than one Nextflow task per orthogroup. Modules therefore take plain `path(files)` collections, and per-file failures are handled inside the shell loop.

**The three-way join** in `workflows/divergence.nf` re-unites each alignment with its tree and manifest by the `<OG>` filename stem (`file.name.tokenize('.')[0]`). `extract_paralogs.py` _always_ writes a `<OG>.tree`, empty when no usable tree exists, so the join is total and no orthogroup is silently dropped. An empty tree file tells `ks.py` to fall back to pairwise-only.

**Failure semantics in CODEML_BATCH** are deliberate: an orthogroup that fails or produces an empty TSV logs a warning and continues; the task only exits non-zero if _every_ orthogroup in the batch failed. `tsv` is `optional: true`. Don't tighten this.

### `bin/ks.py`

The core script. Pairs are formed **within each target species** (`reportable_by_species`), so an ortholog pair can never be emitted — a flat `combinations()` over the manifest would produce them silently. Each row names its `species`.

Produces three independent dS estimates per paralog pair:

- `tree_*` — CODEML M0 (`runmode = 0`) on the whole orthogroup with a fixed topology; dS is the **sum of per-branch dS along the path** between the two paralogs. The primary estimate.
- `pair_*` — CODEML pairwise ML (`runmode = -2`) on just the two sequences. Independent of the tree.
- `yn00_*` — YN00 counting method, a different method class again.

They are not redundant: agreement means the estimate is in a trustworthy regime, and `dS_tree_over_pair` is the cheapest available saturation/alignment/topology alarm. Curation is left to downstream analysis — **no pair is ever filtered out**, quality columns are emitted alongside.

Implementation notes that will bite if forgotten:

- **CODEML is driven through its control file and parsed here**, not via Biopython's `Bio.Phylo.PAML` wrapper. That wrapper's parser raises `IndexError` on truncated result tables, which forced an awkward pairwise fallback in the YN00-based implementation this replaced (removed; see git history for `bin/dnds.py`). Biopython is still used for sequence/alignment/tree handling.
- **Under M0, CODEML prints no `dS tree:` newick.** It prints a per-branch table under `dN & dS for each branch`, with branches labelled by node number (`7..8`). `parse_branch_table` reads it; tips are numbers 1..n in _alignment order_, which is the order `ks.py` writes them, so tips are identified without parsing a tree. Per-branch dS is printed to 4 decimals — a real precision floor for very recent duplicates.
- **Nucleotide matching is exact `record.id` equality** across all CDS files. Strict matching is intentional; a previous attempt at fuzzy ID parsing was reverted because it paired the wrong transcripts. A duplicate ID across two species' CDS files is a hard error, not a silent overwrite.
- **Saturation, not size, is what makes CODEML slow.** A 6-sequence, 319-codon orthogroup of Mycoplasma paralogs (dS ≈ 60) took 44 s for the M0 fit against ~1 s for the same shape of simulated, unsaturated data: the likelihood surface is nearly flat toward large `t` and the optimiser crawls. `params.codeml_timeout` bounds each PAML call (0 = unlimited, the default); a timed-out estimate becomes `NA` and the rest of the orthogroup still runs. The test profile sets `batch_size = 5` and `codeml_timeout = 120` to stay inside its one-hour task cap.
- **Gap stripping is per pair, not per orthogroup.** A single fragmentary paralog would otherwise delete columns from every other pair in the family.
- **`--icode` is the PAML genetic code index**, mapped to NCBI table ids via `PAML_ICODE_TO_NCBI` so the codon alignment is built under the same code CODEML assumes. Getting this wrong silently corrupts both. `params.genetic_code = 3` (Mycoplasma, NCBI table 4) in the test profile.
- Sequences are renamed to short `S00001` ids for PAML's name-length limit and restored afterwards; the PHYLIP is written by hand as strict sequential.
- **The fitted tree is rebuilt from the branch table** by `newick_from_branch_table` and written as `<OG>_dS.nwk` / `<OG>_dN.nwk` with original IDs, including the species not reported on. This is where ortholog distances and internal-node depths live — the TSV only carries target-species paralog pairs. The writer is iterative because a ladder-shaped gene tree would blow the recursion limit.
- **`method` is not `runmode`.** They are separate control-file variables: `runmode` selects the analysis (0 = fit to the supplied tree, −2 = pairwise), `method` selects the branch-length optimiser (0 = all branches simultaneously, 1 = one branch at a time). `params.codeml_method` sets the latter and only affects `runmode = 0`. PAML documents `method = 1` for no-clock analyses, which is what this pipeline always runs; setting `clock = 1` alongside it does not error in this build, so don't rely on codeml to reject the combination. Measured on simulated data: identical dS to 5 d.p. as `method = 0`, ~25× faster at n ≥ 32.

### `bin/extract_paralogs.py`

- **OrthoFinder 3.x writes one combined `Resolved_Gene_Trees/Resolved_Gene_Trees.txt`**, one tree per line as `OG0000000: (newick);` — _not_ a file per orthogroup, and `Gene_Trees/` is empty. `load_combined_trees` handles this; `locate_tree` is the fallback for older per-orthogroup layouts.
- **Gene tree tips are species-prefixed** (`Mycoplasma_agalactiae_gi|290752409|emb|CBH40380.1|`) while `Orthogroups.tsv` and the orthogroup FASTAs use the bare ID. `normalise_tree` strips the prefix, but only when what remains is a gene actually in that orthogroup, so IDs that merely look prefixed are not mangled. `Duplications.tsv` uses the same prefixed form — relevant for downstream LCA work.
- Trees with tips that can't be reconciled are dropped rather than passed on broken.

`bin/generate_cds.py`, `bin/longest_isoform.py`, and `bin/get_matching_cds.py` are standalone data-prep helpers, not called by any process. `generate_cds.py` reverse-translates a proteome into plausible CDS and is how the non-agalactiae test CDS files were made — those are synthetic, so dS values in the test profile are meaningless; the test is a smoke test, not a scientific fixture.

## Containers

`MAFFT_BATCH` and `EXTRACT_PARALOGS` use standard biocontainers. `CODEML_BATCH` uses `docker.io/tdw0student0uml/dnds:1.0`, built from `modules/local/codeml_batch/Dockerfile` (micromamba + biopython 1.84 + paml 4.10.7). **The image already contains `codeml`** — it ships in the same conda `paml` package as `yn00` — so moving off YN00 needed no rebuild. Changing the dependencies does mean rebuilding and pushing, not just editing the conda directive.

## Caveats

- **Two test samplesheets exist** because nf-schema resolves relative paths against the _launch_ directory: `assets/test_samplesheet.csv` (`../../../...`) works under nf-test, `assets/test_samplesheet_manual.csv` (`./...`) works for `nextflow run` from the repo root. `conf/test.config` points at the former, so a direct run needs `--input assets/test_samplesheet_manual.csv`. This is a wart, not a design.
- New scripts in `bin/` **must be `chmod +x`** or Nextflow fails with `bad interpreter: Permission denied`.
- `README.md` and `docs/usage.md` are current. `docs/output.md` and `docs/README.md` are still unmodified nf-core template boilerplate describing a fastq/read-based pipeline with FastQC/MultiQC outputs that do not exist here — don't treat those two as a source of truth.
- `.nf-core.yml` disables a large set of template lint checks (multiqc, igenomes, license, changelog, logos). Don't add the corresponding files back to satisfy a lint error.
- `null/`, `work/`, `diagnose/`, and `tests/out/` are gitignored run detritus. `diagnose/` predates the strict-ID revert and does not reflect current behaviour.
