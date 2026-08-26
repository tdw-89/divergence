```mermaid
flowchart TB
  subgraph DIVERGENCE
    subgraph take
      v0["ch_samplesheet<br/>[meta, fasta, cds]"]
    end

    v1["ch_fastas<br/>(proteomes)"]
    v2["ch_cds<br/>(CDS, collected)"]

    v3([PREPARE_PROTEOME<br/>rename to samplesheet species])
    v4([ORTHOFINDER])
    v5([EXTRACT_PARALOGS])

    v6["*.faa<br/>whole orthogroups"]
    v7["*.tree"]
    v8["*.paralogs.txt"]

    v9([MAFFT_BATCH])
    v10["join on &lt;OG&gt; stem"]
    v11([CODEML_BATCH])

    subgraph emit
      v20["orthofinder"]
      v21["summary"]
      v22["alignments"]
      v23["ks<br/>per-OG + merged ks.tsv"]
    end

    v0 --> v1
    v0 --> v2
    v1 --> v3
    v3 --> v4
    v4 --> v5
    v5 --> v6
    v5 --> v7
    v5 --> v8
    v6 --> v9
    v9 --> v10
    v7 --> v10
    v8 --> v10
    v10 --> v11
    v2 --> v11

    v4 --> v20
    v5 --> v21
    v9 --> v22
    v11 --> v23
  end
```

`PREPARE_PROTEOME` renames each proteome to its samplesheet `species` value: OrthoFinder
takes its species names from filenames, and that name is what `--target_species` matches
and what prefixes gene tree tips. `--orthofinder_dir` skips OrthoFinder (and the rename)
and feeds `EXTRACT_PARALOGS` an existing run directly.

The whole orthogroup goes into the alignment, not just the target species' paralogs:
the extra sequences break up long branches so CODEML's multiple-hit correction has
something to work with. Only pairs listed in `*.paralogs.txt` are reported on.

`EXTRACT_PARALOGS` always writes a `<OG>.tree` — empty when no usable tree exists —
so the three-way join is total and no orthogroup is silently dropped.

The OrthoFinder output directory is published whole, so `Gene_Duplication_Events/`
and `Phylogenetic_Hierarchical_Orthogroups/` are available for last-common-ancestor
and ohnolog analysis downstream. The pipeline itself does not compute those.
