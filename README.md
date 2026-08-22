# divergence

A Nextflow pipeline that estimates pairwise synonymous divergence (dS / Ks) between the paralogs of one or more target species.

## What it does

1. **OrthoFinder** builds orthogroups and gene trees from the input proteomes.
2. **Orthogroup selection** keeps every orthogroup in which at least one target species has at least `--min_paralogs` genes.
3. **MAFFT** aligns each selected orthogroup. The alignment contains the whole orthogroup — all species, not only the targets.
4. **CODEML** (PAML) estimates dS for each pair of paralogs within each target species, using three methods:
   - a maximum-likelihood M0 fit to the whole orthogroup on the OrthoFinder topology, with dS read as the sum of branch dS along the path between the pair;
   - a maximum-likelihood pairwise fit (`runmode = -2`) on the two sequences alone;
   - the Yang & Nielsen (2000) counting method (YN00).

Pairs are formed within a species, so only paralogs are reported, never orthologs. Sequences from non-target species take part in the alignment and the model fit but are not reported on.

The fitted tree is written out with branch lengths in dS and dN units and the original sequence IDs, so any distance the table does not report — between orthologs, or to an internal node — can be read off it.

The pipeline does not filter results. Every pair is emitted with quality columns alongside the estimates.

## Input

A CSV samplesheet with three columns:

```csv title="samplesheet.csv"
species,fasta,cds
Danio_rerio,d_rerio.pep.fa,d_rerio.cds.fa
Oryzias_latipes,o_latipes.pep.fa,o_latipes.cds.fa
Lepisosteus_oculatus,l_oculatus.pep.fa,l_oculatus.cds.fa
```

| Column    | Description                                                                                         |
| --------- | --------------------------------------------------------------------------------------------------- |
| `species` | Species name, no spaces. Target species must be named here exactly as passed to `--target_species`. |
| `fasta`   | Protein FASTA (one sequence per gene).                                                              |
| `cds`     | Nucleotide CDS FASTA. Record IDs must match the protein IDs in `fasta` exactly.                     |

Every species needs a CDS file, not only the target species.

`bin/get_matching_cds.py` produces a CDS file in the required form from a proteome and a genome CDS file. `bin/longest_isoform.py` reduces a proteome to one isoform per gene.

## Usage

```bash
nextflow run . \
   -profile <docker/singularity/conda/...> \
   --input samplesheet.csv \
   --target_species Danio_rerio \
   --outdir <OUTDIR>
```

`--target_species` takes one name, a comma-separated list, or `all`:

```bash
--target_species Danio_rerio,Oryzias_latipes
--target_species all
```

`--input`, `--target_species` and `--outdir` are required. See [`docs/usage.md`](docs/usage.md) for the remaining parameters.

> [!WARNING]
> Provide pipeline parameters on the command line or via `-params-file`. Custom config files supplied with `-c` can set anything **except** parameters.

## Output

```
<OUTDIR>/
├── codeml/       <OG>_ks.tsv        one row per paralog pair
│                 trees/<OG>_dS.nwk  fitted tree, branch lengths in dS
│                 trees/<OG>_dN.nwk  fitted tree, branch lengths in dN
├── mafft/        <OG>.fas           orthogroup alignments
├── extract/      orthogroup_summary.tsv
│                 orthogroups/       per-orthogroup FASTA, tree, paralog manifest
├── orthofinder/  full OrthoFinder output directory
└── pipeline_info/
```

Each row of `<OG>_ks.tsv` is one pair of target-species paralogs:

| Column                                                           | Description                                                       |
| ---------------------------------------------------------------- | ----------------------------------------------------------------- |
| `orthogroup`, `gene_a`, `gene_b`                                 | Pair identity.                                                    |
| `n_seqs_alignment`, `n_codons_alignment`                         | Size of the codon alignment the model was fitted to.              |
| `has_tree`                                                       | Whether a tree-based estimate was possible for this orthogroup.   |
| `n_codons_pair`                                                  | Codons left after dropping columns gapped in either sequence.     |
| `pct_id_a_in_b`, `pct_id_b_in_a`                                 | Nucleotide identity, both directions, gaps counted as mismatches. |
| `tree_dS`, `tree_dN`, `tree_omega`                               | Path sums from the M0 fit.                                        |
| `pair_dS`, `pair_dN`, `pair_omega`, `pair_t`, `pair_S`, `pair_N` | Pairwise CODEML.                                                  |
| `yn00_dS`, `yn00_dN`, `yn00_omega`                               | YN00.                                                             |
| `m0_omega`, `m0_kappa`, `m0_lnL`                                 | Model-level statistics from the M0 fit.                           |
| `dS_tree_over_pair`                                              | `tree_dS / pair_dS`.                                              |

Missing values are `NA`.

`codeml/trees/` holds the tree CODEML fitted for each orthogroup, once with branch lengths in dS units and once in dN. Tips are the original sequence IDs and include the species not reported on, so ortholog distances and internal-node depths can be taken from these directly.

`extract/orthogroup_summary.tsv` lists every orthogroup considered, kept or not, with its sequence and paralog counts, a per-species paralog breakdown, and the reason it was skipped.

The OrthoFinder directory is published in full, including `Gene_Duplication_Events/`, `Phylogenetic_Hierarchical_Orthogroups/` and `Resolved_Gene_Trees/`.

## Credits

divergence was written by Tom Wolfe.

## Citations

Tool references are in [`CITATIONS.md`](CITATIONS.md).

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
