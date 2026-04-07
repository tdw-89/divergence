process MAFFT_BATCH {
    tag "batch"
    label 'process_high'

    conda "bioconda::mafft=7.520"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mafft:7.520--h031d066_3' :
        'quay.io/biocontainers/mafft:7.520--h031d066_3' }"

    input:
    path(fastas)

    output:
    path "*.fas", emit: fas
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    for fasta in ${fastas}; do
        prefix=\$(basename "\$fasta" | sed 's/\\.[^.]*\$//')
        mafft \\
            --thread ${task.cpus} \\
            ${args} \\
            "\$fasta" \\
            > "\${prefix}.fas"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mafft: \$(mafft --version 2>&1 | sed 's/ (.*) //g')
    END_VERSIONS
    """
}
