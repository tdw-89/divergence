process CODEML_BATCH {
    tag "batch"
    label 'process_medium'

    // The published image already contains codeml alongside yn00 (both ship in
    // the paml conda package), so no rebuild was needed to move off YN00.
    conda "bioconda::biopython=1.84 bioconda::paml=4.10.7"
    container 'docker.io/tdw0student0uml/dnds:1.0'

    input:
    path(batch_files)
    path(cds_files)

    output:
    path "*_ks.tsv"    , emit: tsv     , optional: true
    // The fitted tree with branch lengths in dS and dN units, original sequence
    // IDs restored. Carries every distance the TSV does not report, including
    // between orthologs.
    path "*_dS.nwk"    , emit: ds_trees, optional: true
    path "*_dN.nwk"    , emit: dn_trees, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    total=0
    successful=0

    for aln in *.fas; do
        [ -e "\$aln" ] || continue
        og=\${aln%%.*}
        total=\$((total + 1))

        if ks.py \\
            --msa "\$aln" \\
            --cds ${cds_files} \\
            --paralogs "\${og}.paralogs.txt" \\
            --tree "\${og}.tree" \\
            --orthogroup "\$og" \\
            --output "\${og}_ks.tsv" \\
            ${args}; then
            if [[ -s "\${og}_ks.tsv" ]]; then
                successful=\$((successful + 1))
            else
                echo "WARNING: ks.py completed but produced no output for \$og" >&2
            fi
        else
            echo "WARNING: ks.py failed for \$og" >&2
        fi
    done

    # A single orthogroup failing is expected and tolerated; only a batch in
    # which nothing succeeded is treated as a task failure.
    if [[ \${total} -gt 0 && \${successful} -eq 0 ]]; then
        echo "ERROR: CODEML_BATCH produced no dS output in this batch (\${total} orthogroup(s), 0 successful)." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
        paml: \$(codeml 2>&1 | head -1 | sed 's/.*version //;s/,.*//' || echo 'unknown')
    END_VERSIONS
    """
}
