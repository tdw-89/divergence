process DNDS {
    tag "$meta.id"
    label 'process_medium'
    errorStrategy 'ignore'
    
    conda "bioconda::biopython=1.84 bioconda::paml=4.10.7"
    container 'docker.io/tdw0student0uml/dnds:1.0'

    input:
    tuple val(meta), path(msa)
    path nuc_fasta

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    if [[ "${msa}" == *.gz ]]; then
        gunzip -c ${msa} > ${prefix}_msa.fa
        MSA_FILE=${prefix}_msa.fa
    else
        MSA_FILE=${msa}
    fi

    dnds.py \\
        -m \${MSA_FILE} \\
        -f ${nuc_fasta} \\
        -o ${prefix}_dnds.tsv \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        paml: \$(yn00 --version 2>&1 | head -1 || echo 'unknown')
    END_VERSIONS
    """
}
