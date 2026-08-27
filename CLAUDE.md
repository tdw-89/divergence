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

Samplesheet is `species,fasta,cds` (see `assets/test_samplesheet.csv`). **Every species needs a CDS file**, not just the target: codon alignments are built for whole orthogroups. CDS record IDs must match the protein IDs exactly — `bin/get_matching_cds.py` produces files in that form. Input must be uncompressed: OrthoFinder cannot read gzipped FASTA and `ks.py` opens CDS files as text, so `.gz` is rejected by `assets/schema_input.json` rather than failing halfway through.

To re-run the dS analysis without repeating OrthoFinder:

```bash
nextflow run . -profile docker \
  --input samplesheet.csv \
  --orthofinder_dir results/orthofinder/orthofinder \
  --target_species <species_name> \
  --outdir results2
```

Pipeline-level tests (nf-test, `testsDir "."`, profile `test`, nf-core module/subworkflow tests ignored):

```bash
nf-test test                                    # everything
nf-test test tests/default.nf.test              # the pipeline test
nf-test test tests/modules/local/codeml_batch   # the CODEML_BATCH module test
nf-test test --profile test,docker              # override container engine
```

`tests/default.nf.test` reads `orthogroup_summary.tsv` and locates the `status`
column **by name**: it was hardcoded to the wrong index once, and a wrong index counts
zero kept orthogroups rather than failing, which then cascades into the `ks_files.size()`
assertion. Note also that `params` inside an nf-test `then` block resolves against that
test's own `when { params { ... } }` block, **not** `conf/test.config` — anything the
assertions reference (e.g. `target_species`) has to be declared there too.

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

0. **PREPARE_PROTEOME** (local) — copies each proteome to `renamed/<species>.fa`, where `<species>` is the samplesheet's `species` value. **OrthoFinder names species after the filenames it is given**, and that name is the one `--target_species` is resolved against in `Orthogroups.tsv` _and_ the one prefixing every gene tree tip. Without this step the samplesheet's `species` column is read by nf-schema and then never used again, and the real species names are whatever the input files happened to be called — which is how RefSeq accessions full of dots got into tip labels and cost two separate bug fixes. Not published.
1. **ORTHOFINDER** (`modules/nf-core/orthofinder`) — every renamed proteome is `collect()`ed into one channel item tagged `[id: 'orthofinder']`; the prior-run input is always empty. The whole output directory is published, which is what makes downstream LCA/ohnolog analysis possible. **Skipped entirely when `--orthofinder_dir` is given**, which points at a published `orthofinder/` directory from an earlier run: nothing OrthoFinder produces depends on a dS parameter, so re-running the analysis should not repeat the longest step. `--input` stays required either way — the CDS files are read regardless.
2. **EXTRACT_PARALOGS** (local) — runs `bin/extract_paralogs.py`. `--target_species` is a comma-separated list or `all`; each name is resolved against the header by **exact match first, unique prefix second** (errors on zero or multiple matches). Exact-first matters because species now come from the samplesheet, so the header spells them exactly and `human` must not read as ambiguous merely because `human_alt` exists; the prefix fallback survives for `--orthofinder_dir` runs whose species were named from filenames. `resolveTargetSpecies` in `workflows/divergence.nf` applies the same rule **before anything is submitted** — from the samplesheet when we run OrthoFinder, from the supplied directory's header when we do not. That pre-flight check exists because the first real run of this pipeline spent hours on OrthoFinder and then died on a stale `--target_species` that was still spelling species by accession. An orthogroup is kept when any target species has ≥ `--min_paralogs` genes there. For each, emits three files sharing an `<OG>` stem: `<OG>.faa` (the **whole** orthogroup, all species), `<OG>.tree` (gene tree), `<OG>.paralogs.txt` (a `species<TAB>gene` manifest). Also writes `orthogroup_summary.tsv` covering every orthogroup considered, kept or not, including `tree_status` — `not_applicable`, `none_from_orthofinder`, `unparseable`, `unmatched_tips` or `ok`. `has_tree = no` on its own cannot separate "OrthoFinder builds no tree below four sequences" (benign, and most small orthogroups) from "we had a tree and threw it away" (the failure that silently degraded whole runs to pairwise-only, twice).
3. **MAFFT_BATCH** (local) — `.flatten().collate(params.batch_size)` batches the FASTAs; the process loops in bash, one `mafft` call per file, writing `<OG>.fas`. `--anysymbol --localpair --maxiterate 1000` via `ext.args`, i.e. **L-INS-i**: every dS estimate is read off this alignment, so MAFFT's progressive default is not good enough for the divergent paralogs the pipeline exists to measure. `conf/base.config` gives it 6 cores rather than `process_high`'s 12 — L-INS-i wants the memory and wall clock, but thread scaling flattens well before that and only one alignment runs at a time.
4. **CODEML_BATCH** (local) — loops calling `bin/ks.py` once per orthogroup, emitting `<OG>_ks.tsv` plus `<OG>_dS.nwk` / `<OG>_dN.nwk` (published under `codeml/trees/`). The per-orthogroup tables are `collectFile`d into a single `codeml/ks.tsv` for the run; a real run is thousands of files and the concatenation is what downstream analysis wants.

**Why whole orthogroups.** The extra species are never reported on. They exist to break up long branches: a paralog pair whose duplication predates several speciations sits at the end of one enormous branch, where dS is saturated and the multiple-hit correction has nothing to work with. Interposing orthologs subdivides that path into individually short branches. Do not "optimise" this back to target-species-only alignments — it is the point of the design.

**Batching** — orthogroup counts run into the thousands, so `params.batch_size` (default 20) controls how many orthogroups each MAFFT/CODEML task handles in a bash loop, rather than one Nextflow task per orthogroup. Modules therefore take plain `path(files)` collections, and per-file failures are handled inside the shell loop.

**Batches are packed by cost, not by channel order** (`balancedBatches` in `workflows/divergence.nf`). This matters more than it looks: channel order is filename order, and **OrthoFinder numbers orthogroups by descending size**, so a plain `.collate()` hands the first task the twenty largest orthogroups in the entire run. Measured on a real primate run that is **40% of all pairwise work in one task**, with a max/median task cost ratio of 761x — one task ran past an eight-hour wall clock and was killed (SLURM exit 140) while the median finished in 28 s. Dealing a cost-sorted list round-robin brings the worst task to ~4% and the ratio to 12x. The cost model (`manifestPairCost`) is the exact pair count from the manifest, capped at `max_pairwise_pairs`, plus `n^1.6 / 8` for the tree fit; it only has to get the _ordering_ roughly right. MAFFT batches use file size as the proxy.

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
- **CODEML's `MAXNSONS` is 3, and OrthoFinder emits wider nodes freely.** A node with four or more daughters makes codeml abort with `error: too many daughter nodes, raise MAXNSONS` in under a second, which costs the whole orthogroup its `tree_dS` — the primary estimate. `_resolve_polytomies` splits multifurcations into binary nodes joined by zero-length branches before the tree is written. Polytomies are an artefact of collapsed zero-length branches in the first place, so re-imposing an arbitrary order is as defensible as the collapse was, and CODEML re-estimates every branch anyway. Resolution is **balanced, not a ladder**: a ladder is as deep as the polytomy is wide, and Biopython walks trees recursively, so a few thousand paralogs under one node would overflow the stack. Measured on a real primate run: 448 of 2,354 orthogroups (19%) lost their tree to this, every failure having max degree ≥ 4 and every success exactly 3.
- **The codon alignment drops bad sequences rather than failing the orthogroup.** `Bio.codonalign.build` is all or nothing: one protein its CDS does not translate back to — an annotated frameshift, a partial CDS, a selenocysteine — raises and takes the entire orthogroup with it. RefSeq carries these at a percent or two of sequences, and at `1-(1-p)^n` that destroyed 26% of kept orthogroups in a real primate run, rising monotonically with size (7% at 2–4 sequences, 66% at 65+) — i.e. it deleted exactly the large families. `build_codon_alignment_tolerantly` tries the whole orthogroup first (so a clean one costs nothing extra) and only then screens sequence by sequence, using Biopython's own `build` on a single record as the test rather than reimplementing its notion of a match. Dropped ids are warned about, and `retained`/`reportable_by_species`/the gene tree are all re-derived afterwards — they must agree on which sequences survived.
- **`run_codeml` returns a `PamlRun`, not a string.** PAML reports fatal errors on stdout and _still_ leaves a partial `mlc` behind, so reading only that file gives a truncated result with no reason attached. Both bugs above were invisible for exactly this reason. `PamlRun.diagnostic()` pulls the `error:` line out of the console log and every parse failure quotes it.
- **Nucleotide matching is exact `record.id` equality** across all CDS files. Strict matching is intentional; a previous attempt at fuzzy ID parsing was reverted because it paired the wrong transcripts. A duplicate ID across two species' CDS files is a hard error, not a silent overwrite. The lookup is `CdsIndex`, a `Mapping` over one `SeqIO.index` per file: an orthogroup needs a few dozen of the run's hundreds of thousands of records, and offsets cost a fraction of parsed `SeqRecord`s. It holds file handles, so `main` closes it once the codon alignment is built.
- **Saturation, not size, is what makes CODEML slow.** A 6-sequence, 319-codon orthogroup of Mycoplasma paralogs (dS ≈ 60) took 44 s for the M0 fit against ~1 s for the same shape of simulated, unsaturated data: the likelihood surface is nearly flat toward large `t` and the optimiser crawls. `params.codeml_timeout` bounds each PAML call (0 = unlimited, the default); a timed-out estimate becomes `NA` and the rest of the orthogroup still runs. The test profile sets `batch_size = 5` and `codeml_timeout = 120` to stay inside its one-hour task cap.
- **The cross-checks are sampled above `params.max_pairwise_pairs` (default 5000).** `tree_dS` comes from one M0 fit per orthogroup and covers every pair by path summation, so the primary estimate's cost grows with orthogroup _size_. `pair_*` and `yn00_*` cost two external processes per pair — measured at **65 ms/pair and near-constant in alignment width** (200 to 3000 codons all cost the same, so it is process spawn and file I/O, not likelihood evaluation) — and pairs grow quadratically. A 417-paralog family is 86,751 pairs, ~1.5 h of spawning, to check a fit that took minutes. Sampling is seeded on the orthogroup name so reruns pick the same pairs, and `crosschecked` distinguishes a sampled-out `NA` from a failed one. **`crosscheck_cap_applies` keys the cap on orthogroup size, not on whether the tree fit worked.** Keying it on success — as it did originally — meant a timeout, a `MAXNSONS` abort or an unreconcilable tip _lifted_ the cap, and those happen to the largest families, so the orthogroups least able to afford the full quadratic cross-check were the only ones that ran it. `manifestPairCost` compounds it by predicting the capped cost for them, so `balancedBatches` underestimates exactly those tasks by up to the cap ratio. The cap is lifted only below `MIN_TREE_SEQS` (4) sequences, where OrthoFinder built no tree in the first place and there are at most three pairs.
- **`params.codeml_threads` splits each CODEML_BATCH task's cores.** The module runs `task.cpus / codeml_threads` orthogroups at a time under `xargs -P`, each `ks.py` running that many cross-checks in a `ThreadPoolExecutor`. The two never oversubscribe. Threads are enough because the work is in external processes. Default 1: raise it only when batches are dominated by a few very large families, where there are too few orthogroups to fill the cores. Note the interaction with `codeml_timeout`, which is wall-clock — fewer concurrent workers means each PAML call gets more CPU and is _less_ likely to time out, which shows up as slower runs that produce more results.
- **Pair coverage is `n_codons_pair / min(n_codons_a, n_codons_b)`, never over `n_codons_alignment`.** The alignment width is the union of every indel in the family, so dividing by it makes a healthy orthogroup of variable-length proteins look identical to a broken one — this metric produced a false diagnosis of "MCL-chained orthogroups" that a re-measure against the shorter sequence demolished: across 411 primate orthogroups including every large one, the median orthogroup scored **1.00** and **none** fell below 0.5. The real pathology is per-_pair_ fragments and partial annotations (a 364-codon sequence sharing literally zero columns with an 896-codon one), which is why the per-sequence lengths are emitted. **It is not partitionable**: the "shares alignable data" relation is non-transitive (measured 136/502 violations in one orthogroup), so no subdivision of the orthogroup — no tree cut, no re-clustering — can express it. Filter on the column instead.
- **Gap stripping is per pair, not per orthogroup.** A single fragmentary paralog would otherwise delete columns from every other pair in the family.
- **`--icode` is the PAML genetic code index**, mapped to NCBI table ids via `PAML_ICODE_TO_NCBI` so the codon alignment is built under the same code CODEML assumes. Getting this wrong silently corrupts both. `params.genetic_code = 3` (Mycoplasma, NCBI table 4) in the test profile.
- Sequences are renamed to short `S00001` ids for PAML's name-length limit and restored afterwards; the PHYLIP is written by hand as strict sequential.
- **The fitted tree is rebuilt from the branch table** by `newick_from_branch_table` and written as `<OG>_dS.nwk` / `<OG>_dN.nwk` with original IDs, including the species not reported on. This is where ortholog distances and internal-node depths live — the TSV only carries target-species paralog pairs. The writer is iterative because a ladder-shaped gene tree would blow the recursion limit.
- **`method` is not `runmode`.** They are separate control-file variables: `runmode` selects the analysis (0 = fit to the supplied tree, −2 = pairwise), `method` selects the branch-length optimiser (0 = all branches simultaneously, 1 = one branch at a time). `params.codeml_method` sets the latter and only affects `runmode = 0`. PAML documents `method = 1` for no-clock analyses, which is what this pipeline always runs; setting `clock = 1` alongside it does not error in this build, so don't rely on codeml to reject the combination. Measured on simulated data: identical dS to 5 d.p. as `method = 0`, ~25× faster at n ≥ 32.

#### Measured cost

From this repo's own benchmarks, so the batching and cap parameters can be set from data rather than guessed. Re-measure rather than trusting these if the container's PAML version changes.

- **M0 tree fit** (`method = 1`, 300 codons, dS ≈ 1, unrooted, one replicate per point): n = 64 → 18.6 s, 128 → 50.0 s, 256 → 259 s, 512 → 530 s. That is 28.5× for 8× the sequences, so **≈ O(n^1.6)** over that range. The per-doubling exponents bounce (1.43, 2.38, 1.03) because codeml's iteration count depends on the particular data realisation — read the endpoints as signal, the segments as noise.
- **Pairwise + YN00 end to end through `ks.py`** (500 codons): **0.088 s/pair**, linear in pair count and flat in divergence — raw `runmode = -2` is 0.07 s at dS ≈ 0.5 and 0.09 s at dS ≈ 4. Measured as 190 pairs / 17 s and 780 pairs / 68 s.
- **So paralog count is cheap and orthogroup size is not.** `C(100,2) = 4950` pairs is ~7 min; a 500-sequence tree fit is ~9 min at 300 codons, ~15 min at vertebrate CDS length. `params.codeml_timeout` bounds the tree fit but is no defence against pair count: each pairwise call is already sub-0.1 s, so only the count grows, and the count is cheap. `params.max_pairwise_pairs` is what bounds it.
- **Re-indexing every CDS file on every `ks.py` call is not a bottleneck**: 176,000 records across 8 files (297 MB) parse in 1.9 s at 0.54 GB RSS. That was the measurement behind moving to `SeqIO.index` anyway — the time was tolerable but the half-gigabyte was paid again by every concurrent `xargs` worker.
- **Caveat**: these simulations mutate synonymous sites only (ω ≈ 0), the easy case for the optimiser. Real ω heterogeneity, alignment error and true saturation are all slower — see the Mycoplasma dS ≈ 60 case above.

### `bin/extract_paralogs.py`

- **OrthoFinder 3.x writes one combined `Resolved_Gene_Trees/Resolved_Gene_Trees.txt`**, one tree per line as `OG0000000: (newick);` — _not_ a file per orthogroup, and `Gene_Trees/` is empty. `load_combined_trees` handles this; `locate_tree` is the fallback for older per-orthogroup layouts.
- **OrthoFinder spells the species prefix differently in the two files it writes.** `Orthogroups.tsv` keeps the name verbatim; gene-tree tips have _some_ punctuation rewritten to underscores. Which characters is not documented and not worth guessing — observed output rewrites `.` but leaves `-` alone, so `GCF_053564925.1_Olatipes_Hd-rR_3.1` becomes `GCF_053564925_1_Olatipes_Hd-rR_3_1` on a tip. `prefix_pattern` therefore builds a regex in which every non-alphanumeric character of the species name matches either itself or `_`, rather than hard-coding a substitution. Over-matching is safe because `normalise_tree` only strips a prefix when the remainder is a gene genuinely in that orthogroup. Getting this wrong drops **every** tree for a run whose species names contain punctuation (RefSeq accessions do), silently degrading every dS estimate to pairwise-only — `has_tree` is `no` across the board and `tree_dS` is `NA`. Two separate bugs here have had exactly that effect: first matching only the verbatim spelling, then over-correcting by rewriting every non-alphanumeric character, which missed the one species with a hyphen. **One unmatched tip drops the whole tree**, so a single mis-spelt species takes down every orthogroup it appears in — which is most of them. The gene id after the prefix is never rewritten, so only the species side is treated this way. Since `PREPARE_PROTEOME` now names species from the samplesheet, a run can avoid the whole problem by using only `[A-Za-z0-9_]` there — but `prefix_pattern` stays, because nothing enforces that and `--orthofinder_dir` can supply a run built any way at all.
- **Gene tree tips are species-prefixed** (`Mycoplasma_agalactiae_gi|290752409|emb|CBH40380.1|`) while `Orthogroups.tsv` and the orthogroup FASTAs use the bare ID. `normalise_tree` strips the prefix, but only when what remains is a gene actually in that orthogroup, so IDs that merely look prefixed are not mangled. `Duplications.tsv` uses the same prefixed form — relevant for downstream LCA work.
- Trees with tips that can't be reconciled are dropped rather than passed on broken.

`bin/generate_cds.py`, `bin/longest_isoform.py` and `bin/get_matching_cds.py` are standalone data-prep helpers, not called by any process. `longest_isoform.py` has two schemes because annotation sources put the gene identifier in different places: by default it regexes one out of the protein defline (Ensembl-style), while **RefSeq protein deflines carry no gene id at all**, so `--refseq` drives selection from the CDS file instead, grouping on `[db_xref=GeneID:...]`, keeping the longest CDS per gene and pulling its protein out by `protein_id`. `--all` applies either scheme across a directory of per-species subdirectories, naming outputs `<input>_longest.<ext>` and tolerating a per-species failure. `generate_cds.py` reverse-translates a proteome into plausible CDS and is how the non-agalactiae test CDS files were made — those are synthetic, so dS values in the test profile are meaningless; the test is a smoke test, not a scientific fixture.

## Containers

`MAFFT_BATCH` and `EXTRACT_PARALOGS` use standard biocontainers. `CODEML_BATCH` uses `docker.io/tdw0student0uml/dnds:1.0`, built from `modules/local/codeml_batch/Dockerfile` (micromamba + biopython 1.84 + paml 4.10.7). **The image already contains `codeml`** — it ships in the same conda `paml` package as `yn00` — so moving off YN00 needed no rebuild. Changing the dependencies does mean rebuilding and pushing, not just editing the conda directive.

## Caveats

- **Two test samplesheets exist** because nf-schema resolves relative paths against the _launch_ directory: `assets/test_samplesheet.csv` (`../../../...`) works under nf-test, `assets/test_samplesheet_manual.csv` (`./...`) works for `nextflow run` from the repo root. `conf/test.config` points at the former, so a direct run needs `--input assets/test_samplesheet_manual.csv`. This is a wart, not a design.
- New scripts in `bin/` **must be `chmod +x`** or Nextflow fails with `bad interpreter: Permission denied`.
- `README.md` and `docs/usage.md` are current. `docs/output.md` and `docs/README.md` are still unmodified nf-core template boilerplate describing a fastq/read-based pipeline with FastQC/MultiQC outputs that do not exist here — don't treat those two as a source of truth.
- `.nf-core.yml` disables a large set of template lint checks (multiqc, igenomes, license, changelog, logos). Don't add the corresponding files back to satisfy a lint error.
- **`params.codeml_timeout` is wall-clock, and `CODEML_BATCH` now runs `xargs -P ${task.cpus}`**, so contention between concurrent workers eats into each call's budget. Observed once locally (10 timeouts at `codeml_timeout = 120` under 4-way concurrency on emulated amd64, against 0 in an earlier serial run) but **not confirmed against a controlled serial baseline** — the mechanism is plausible, the magnitude is not established. It matters most where cores are oversubscribed: `conf/test.config` pairs `resourceLimits cpus: 4` with `codeml_timeout = 120`, so a small CI runner may time out more than a dev box. Timeouts degrade gracefully (`tree_dS` becomes `NA`, pairwise and YN00 still populate), so the symptom is thinner tree coverage rather than failure.
- **`MAFFT_BATCH` has no per-file failure tolerance and no timeout**, unlike `CODEML_BATCH`: the loop runs under `bash -ue`, so one failed `mafft` fails the whole batch. `params.max_orthogroup_seqs` is applied in `extract_paralogs.py` before anything is written, which makes it the only guard on the MAFFT side — and it matters more now that MAFFT runs L-INS-i, which is quadratic in sequence count.
- `null/`, `work/`, `diagnose/`, and `tests/out/` are gitignored run detritus. `diagnose/` predates the strict-ID revert and does not reflect current behaviour.
