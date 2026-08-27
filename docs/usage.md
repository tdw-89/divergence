# divergence: Usage

## Samplesheet input

`--input` takes a CSV with a header row and three columns:

```csv title="samplesheet.csv"
species,fasta,cds
Danio_rerio,d_rerio.pep.fa,d_rerio.cds.fa
Oryzias_latipes,o_latipes.pep.fa,o_latipes.cds.fa
Lepisosteus_oculatus,l_oculatus.pep.fa,l_oculatus.cds.fa
```

| Column    | Requirement                                                         |
| --------- | ------------------------------------------------------------------- |
| `species` | No spaces. The target species' value must match `--target_species`. |
| `fasta`   | Protein FASTA, extension `.fa`, `.faa` or `.fasta`.                 |
| `cds`     | Nucleotide CDS FASTA, extension `.fa`, `.fna` or `.fasta`.          |

Requirements on the sequences themselves:

- One row per species. Every species needs both files.
- Uncompressed only. OrthoFinder cannot read gzipped FASTA.
- The `species` value names the species everywhere downstream: the proteome is renamed to it before OrthoFinder runs, so it is what `--target_species` matches and what prefixes gene tree tips. Stick to letters, digits and underscores — OrthoFinder rewrites some punctuation in tip labels and not in the orthogroup table.
- Every record ID in `cds` must match a protein ID in the corresponding `fasta` **exactly**. No version-suffix or substring matching is performed. Protein IDs with no CDS match are dropped from the analysis with a warning.
- Record IDs must be unique across all CDS files. A duplicate ID in two species' files is a fatal error.
- One isoform per gene. Multiple isoforms of the same gene will be treated as paralogs.
- CDS must be in frame and translate to the corresponding protein under the genetic code given by `--genetic_code`.

Relative paths in the samplesheet are resolved against the directory you launch from, not the samplesheet's location.

### Preparing the input files

Two helper scripts in `bin/` are not part of the pipeline but produce input in the required form. Outputs are written beside their inputs as `<input>_longest.<ext>` unless a path is given:

```bash
# reduce a proteome to the longest isoform per gene (gene id read from the
# protein defline, e.g. Ensembl's `gene:ENSDARG...`)
longest_isoform.py --input proteome.fa

# NCBI RefSeq: protein deflines carry no gene id, so selection is driven from
# the CDS file and both files are written with matching record IDs
longest_isoform.py --refseq --input GCF_x_protein.faa --cds GCF_x_cds_from_genomic.fna

# ...or across a directory holding one subdirectory per species
longest_isoform.py --refseq --all ./genomes

# pull the CDS matching each protein; a CDS is accepted only if the protein ID
# appears in its defline and it translates exactly to that protein
get_matching_cds.py proteome_longest.fa genome_cds.fa species_cds.fa
```

`longest_isoform.py` finds the gene ID with the regex given by `-r`, which defaults to
`gene:\S+` and so expects Ensembl-style deflines. Pass a different pattern for other
sources.

## Running the pipeline

```bash
nextflow run . \
   -profile docker \
   --input samplesheet.csv \
   --target_species Danio_rerio \
   --outdir results
```

To analyse several species in one run, pass a comma-separated list, or `all` for every
species in the samplesheet:

```bash
--target_species Danio_rerio,Oryzias_latipes
--target_species all
```

`--target_species` takes the samplesheet's `species` values, not filenames or accessions:
the proteomes are renamed to those values before OrthoFinder runs, so they are what the
orthogroup table is keyed on. The names are checked before anything is submitted, so a
stale value fails in seconds instead of after OrthoFinder. (A name that is a prefix of
another species' name still resolves to itself; only a name matching several species and
none of them exactly is an error.)

Pairs are always formed within a species, so adding target species never produces
ortholog pairs. It costs only the extra pairwise fits: the orthogroup alignment and the
tree fit are shared across all target species in that orthogroup.

To check the setup before a real run:

```bash
nextflow run . -profile test,docker --outdir results
```

Resume an interrupted run with `-resume`.

### Reusing an OrthoFinder run

OrthoFinder is by far the longest step, and nothing it produces depends on the dS
parameters. To re-run the dS analysis without repeating it, point `--orthofinder_dir` at
the `orthofinder/` directory an earlier run published:

```bash
nextflow run . -profile docker \
   --input samplesheet.csv \
   --orthofinder_dir results/orthofinder/orthofinder \
   --target_species Danio_rerio \
   --outdir results2
```

`--input` is still required — the CDS files are read either way — and the samplesheet's
species names must be the ones that run was built with.

### Required parameters

| Parameter          | Description                                                                                                                                                                                                                                   |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input`          | Path to the samplesheet.                                                                                                                                                                                                                      |
| `--target_species` | Species whose paralogs are analysed: one name, a comma-separated list, or `all`. Each must be a `species` value from the samplesheet. Resolved before any work is submitted, so a wrong name fails immediately rather than after OrthoFinder. |
| `--outdir`         | Output directory.                                                                                                                                                                                                                             |

### Optional inputs

| Parameter           | Description                                                                                |
| ------------------- | ------------------------------------------------------------------------------------------ |
| `--orthofinder_dir` | Reuse the `orthofinder/` directory of an earlier run instead of running OrthoFinder again. |

### Analysis parameters

| Parameter               | Default | Description                                                                                                                                   |
| ----------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `--genetic_code`        | `0`     | PAML genetic code index (`icode`). `0` universal, `3` Mycoplasma/Spiroplasma (NCBI table 4). Used for both the codon alignment and CODEML.    |
| `--min_paralogs`        | `2`     | Minimum target-species genes for an orthogroup to be analysed.                                                                                |
| `--max_orthogroup_seqs` | `0`     | Skip orthogroups larger than this in total sequences. `0` disables.                                                                           |
| `--batch_size`          | `20`    | Orthogroups handled per MAFFT or CODEML task.                                                                                                 |
| `--codeml_method`       | `1`     | CODEML `method`. `0` optimises all branches at once, `1` one at a time. `1` is much faster on large orthogroups and gives the same estimates. |
| `--codeml_timeout`      | `0`     | Abandon a single CODEML or YN00 call after this many seconds. `0` waits indefinitely.                                                         |
| `--max_pairwise_pairs`  | `5000`  | Cross-check at most this many pairs per orthogroup, sampled at random. `0` checks all.                                                        |
| `--codeml_threads`      | `1`     | Pairwise cross-checks run concurrently within one orthogroup.                                                                                 |
| `--skip_pairwise`       | `false` | Omit the pairwise CODEML estimate.                                                                                                            |
| `--skip_yn00`           | `false` | Omit the YN00 estimate.                                                                                                                       |

### Judging a pair

Coverage is `n_codons_pair / min(n_codons_a, n_codons_b)`. Do not divide by
`n_codons_alignment`: that is the union of every indel in the family, so a perfectly
healthy orthogroup of variable-length proteins scores as badly as one containing a
fragment. Measured against the shorter sequence, real orthogroups sit at ~1.0 (across a
real primate run, the median orthogroup scored exactly 1.00 and none fell below 0.5),
while individual pairs in which one sequence is a partial annotation drop to near 0 --
sometimes literally zero shared codons between a 364-codon fragment and a 896-codon gene.

That is a property of a _pair_, not of a subfamily: the same orthogroup usually contains
pairs that overlap perfectly. It cannot be fixed by splitting the orthogroup, and it does
not need to be -- filter on the column.

### Cost, and where it actually goes

The tree-based `tree_dS` is the primary estimate, and it comes from **one** M0 fit per
orthogroup that yields a value for every pair by summing along paths. Its cost grows with
orthogroup size, not with paralog count — a 700-tip fit is a few minutes.

The `pair_*` and `yn00_*` columns are independent cross-checks, and they cost two external
processes **per pair** — about 65 ms, near-constant regardless of alignment length. Pairs
grow quadratically: a 417-paralog family is 86,751 pairs, roughly an hour and a half of
process spawning, to check a fit that took minutes. `--max_pairwise_pairs` caps that by
sampling. Rows outside the sample keep their `tree_dS` and report `crosschecked = no`, so
an `NA` from sampling is distinguishable from one caused by a failure. The cap is lifted
only for orthogroups too small for OrthoFinder to have built a gene tree (under four
sequences), where the pairwise columns are the only result and there are at most three
pairs; an orthogroup whose tree fit failed stays capped, because those are the large
families. Sampling is seeded on the orthogroup name, so reruns pick the same pairs.

MAFFT runs L-INS-i (`--localpair --maxiterate 1000`), which is quadratic in sequence
count: every dS estimate is read off that alignment, so the progressive default is not
accurate enough for divergent paralogs. `--max_orthogroup_seqs` is the only guard on that
side — MAFFT has no per-file timeout and one failure fails its batch.

Sequences diverged past saturation make CODEML's optimiser crawl after an unbounded dS,
which can make the tree fit slow regardless of orthogroup size. `--codeml_timeout` bounds
it: a timed-out estimate is reported as `NA`, a warning is logged, and the rest of the
orthogroup still runs. Note that the timeout is wall-clock, so concurrent workers within a
task compete for it.

Batches are packed by estimated cost rather than in orthogroup order. OrthoFinder numbers
orthogroups by descending size, so naive batching would hand the first task the largest
twenty in the run — measured at 40% of a real dataset's pairwise work in one task. Cost-
balanced packing brings that to about 4%.

`--codeml_threads` divides each task's cores between orthogroups and pairs: the task runs
`task.cpus / codeml_threads` orthogroups at a time, each using that many threads, so the
two never oversubscribe. Leave it at `1` unless batches are dominated by a few very large
families, where there are too few orthogroups to spread across the cores.

Full parameter documentation is generated from the schema:

```bash
nextflow run . --help
nextflow run . --help_full
```

## Profiles

Use `-profile` to select a container engine. More than one profile can be given, comma-separated and in order of precedence.

| Profile                                                                   | Effect                                    |
| ------------------------------------------------------------------------- | ----------------------------------------- |
| `test`                                                                    | Small bundled dataset.                    |
| `docker`, `singularity`, `apptainer`, `podman`, `shifter`, `charliecloud` | Container engine.                         |
| `conda`, `mamba`                                                          | Conda environments instead of containers. |
| `emulate_amd64`                                                           | Run amd64 images on an arm64 host.        |
| `unity`                                                                   | UMass Unity cluster configuration.        |
| `wave`, `arm64`, `gpu`, `debug`                                           | See `nextflow.config`.                    |

## Resources and configuration

Process resources come from labels in `conf/base.config`. To change them, supply a custom config with `-c`:

```groovy title="custom.config"
process {
    withName: 'CODEML_BATCH' {
        cpus   = 8
        memory = 32.GB
        time   = 24.h
    }
}
```

```bash
nextflow run . -profile docker -c custom.config --input samplesheet.csv \
  --target_species Danio_rerio --outdir results
```

Parameters cannot be set from a `-c` config. Use the command line or `-params-file`.

## Output beyond the tables

`codeml/ks.tsv` is every pair in the run in one table — the concatenation of the
per-orthogroup `codeml/<OG>_ks.tsv` files, which are also kept.

`codeml/trees/<OG>_dS.nwk` and `<OG>_dN.nwk` hold the tree CODEML fitted for each
orthogroup, with branch lengths in dS and dN units and the original sequence IDs on the
tips. They include the species not reported on, so distances the tables omit — between
orthologs, or from a tip to an internal node — can be read off them. The `tree_dS` column
of the table is the sum of dS branch lengths along the path between the two paralogs in
this tree.

## Failure behaviour

An orthogroup whose dS estimation fails logs a warning and is skipped; the run continues. A CODEML task fails only when every orthogroup in its batch failed. Orthogroups with no usable gene tree still produce pairwise and YN00 estimates, with the `tree_*` columns set to `NA`.

## Requirements

- Nextflow `>=25.04.0`, running on Java 17–24.
- A container engine, or conda.
