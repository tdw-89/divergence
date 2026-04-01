process EXTRACT_PARALOGS {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::biopython=1.84"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.84' :
        'biocontainers/biopython:1.84' }"

    input:
    tuple val(meta), path(ortho_dir)
    val species_name

    output:
    // We output all generated FASTA files. 
    // They don't have a meta map yet, we will add that in the workflow.
    path "*_unaligned.fasta", emit: unaligned_fastas
    path "versions.yml"     , emit: versions

    script:
    """
    extract_paralogs.py \\
        --tsv ${ortho_dir}/Orthogroups/Orthogroups.tsv \\
        --seqdir ${ortho_dir}/Orthogroup_Sequences \\
        --species ${species_name}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}