process DNDS_BATCH {
    tag "batch"
    label 'process_medium'

    conda "bioconda::biopython=1.84 bioconda::paml=4.10.7"
    container 'docker.io/tdw0student0uml/dnds:1.0'

    input:
    path(msas)
    path nuc_fasta

    output:
    path "*_dnds.tsv", emit: tsv, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    total_msas=0
    successful_msas=0

    for msa in ${msas}; do
        ((total_msas+=1))
        prefix=\$(basename "\$msa" | sed 's/\\.[^.]*\$//')

        if [[ "\$msa" == *.gz ]]; then
            gunzip -c "\$msa" > "\${prefix}_decompressed.fa"
            actual_msa="\${prefix}_decompressed.fa"
        else
            actual_msa="\$msa"
        fi

        if dnds.py \\
            -m "\${actual_msa}" \\
            -f ${nuc_fasta} \\
            -o "\${prefix}_dnds.tsv" \\
            ${args}; then
            if [[ -s "\${prefix}_dnds.tsv" ]]; then
                ((successful_msas+=1))
            else
                echo "WARNING: dnds.py completed but produced no output for \$msa" >&2
            fi
        else
            echo "WARNING: dnds.py failed for \$msa" >&2
        fi
    done

    if [[ \${total_msas} -gt 0 && \${successful_msas} -eq 0 ]]; then
        echo "ERROR: DNDS_BATCH produced no dN/dS outputs in this batch (\${total_msas} input MSA(s), 0 successful)." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        paml: \$(yn00 --version 2>&1 | head -1 || echo 'unknown')
    END_VERSIONS
    """
}