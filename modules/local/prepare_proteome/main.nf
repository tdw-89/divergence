process PREPARE_PROTEOME {
    tag "$meta.id"
    label 'process_single'

    conda "conda-forge::coreutils=9.5"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:22.04' :
        'nf-core/ubuntu:22.04' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("renamed/${meta.id}.fa"), emit: fasta

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # OrthoFinder names each species after the *filename* it was given, and that
    # name is what every downstream step matches on: the Orthogroups.tsv header
    # that --target_species is resolved against, and the species prefix on gene
    # tree tips. Renaming the proteome to the samplesheet's species name here is
    # what makes the samplesheet the single source of truth for species naming.
    #
    # Written into a subdirectory so the output never collides with the staged
    # input when a file is already correctly named.
    mkdir renamed
    cp -L ${fasta} renamed/${meta.id}.fa
    """

    stub:
    """
    mkdir renamed
    touch renamed/${meta.id}.fa
    """
}
