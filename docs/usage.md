# divergence: Usage

## Samplesheet input

`--input` takes a CSV with a header row and three columns:

```csv title="samplesheet.csv"
species,fasta,cds
Danio_rerio,d_rerio.pep.fa,d_rerio.cds.fa
Oryzias_latipes,o_latipes.pep.fa,o_latipes.cds.fa
Lepisosteus_oculatus,l_oculatus.pep.fa,l_oculatus.cds.fa
```

| Column    | Requirement                                                                 |
| --------- | --------------------------------------------------------------------------- |
| `species` | No spaces. The target species' value must match `--target_species`.         |
| `fasta`   | Protein FASTA, extension `.fa`, `.faa`, `.fasta` (optionally `.gz`).        |
| `cds`     | Nucleotide CDS FASTA, extension `.fa`, `.fna`, `.fasta` (optionally `.gz`). |

Requirements on the sequences themselves:

- One row per species. Every species needs both files.
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

Pairs are always formed within a species, so adding target species never produces
ortholog pairs. It costs only the extra pairwise fits: the orthogroup alignment and the
tree fit are shared across all target species in that orthogroup.

To check the setup before a real run:

```bash
nextflow run . -profile test,docker --outdir results
```

Resume an interrupted run with `-resume`.

### Required parameters

| Parameter          | Description                                                                                                                                                                                           |
| ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--input`          | Path to the samplesheet.                                                                                                                                                                              |
| `--target_species` | Species whose paralogs are analysed: one name, a comma-separated list, or `all`. Each name is matched as a prefix against the OrthoFinder orthogroup table header and must match exactly one species. |
| `--outdir`         | Output directory.                                                                                                                                                                                     |

### Analysis parameters

| Parameter               | Default | Description                                                                                                                                   |
| ----------------------- | ------- | --------------------------------------------------------------------------------------------------------------------------------------------- |
| `--genetic_code`        | `0`     | PAML genetic code index (`icode`). `0` universal, `3` Mycoplasma/Spiroplasma (NCBI table 4). Used for both the codon alignment and CODEML.    |
| `--min_paralogs`        | `2`     | Minimum target-species genes for an orthogroup to be analysed.                                                                                |
| `--max_orthogroup_seqs` | `0`     | Skip orthogroups larger than this in total sequences. `0` disables.                                                                           |
| `--max_target_paralogs` | `0`     | Skip orthogroups with more target-species genes than this. `0` disables.                                                                      |
| `--batch_size`          | `20`    | Orthogroups handled per MAFFT or CODEML task.                                                                                                 |
| `--codeml_method`       | `0`     | CODEML `method`. `0` optimises all branches at once, `1` one at a time. `1` is much faster on large orthogroups and gives the same estimates. |
| `--codeml_timeout`      | `0`     | Abandon a single CODEML or YN00 call after this many seconds. `0` waits indefinitely.                                                         |
| `--skip_pairwise`       | `false` | Omit the pairwise CODEML estimate.                                                                                                            |
| `--skip_yn00`           | `false` | Omit the YN00 estimate.                                                                                                                       |

Pairwise cost grows with the square of the paralog count, so a single large gene family
can dominate a run; `--max_target_paralogs` caps it. Separately, sequences diverged past
saturation make CODEML's optimiser crawl after an unbounded dS, which can make the tree
fit slow regardless of orthogroup size. `--codeml_timeout` bounds it: a timed-out estimate
is reported as `NA`, a warning is logged, and the rest of the orthogroup still runs.

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
