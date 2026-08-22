process EXTRACT_PARALOGS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::biopython=1.84"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(ortho_dir)
    val species_names

    output:
    // Whole orthogroups, not just the target species' genes: the other species
    // break up long branches so the codon model can correct for multiple hits.
    // These channels are bare paths, keyed to each other by the <OG> filename
    // stem, and are re-joined in the workflow.
    path "orthogroups/*.faa"          , emit: fastas
    path "orthogroups/*.tree"         , emit: trees
    path "orthogroups/*.paralogs.txt" , emit: manifests
    path "orthogroup_summary.tsv"     , emit: summary
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    extract_paralogs.py \\
        --tsv ${ortho_dir}/Orthogroups/Orthogroups.tsv \\
        --seqdir ${ortho_dir}/Orthogroup_Sequences \\
        --treedir ${ortho_dir}/Resolved_Gene_Trees \\
        --species '${species_names}' \\
        --outdir orthogroups \\
        ${args}

    mv orthogroups/orthogroup_summary.tsv orthogroup_summary.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}
