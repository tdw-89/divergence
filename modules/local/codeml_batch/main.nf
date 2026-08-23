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
    # ks.py is single-threaded and PAML has no threading of its own, so the
    # orthogroups in a batch are run concurrently rather than in sequence.
    # Each ks.py call already makes its own tempfile.mkdtemp() and runs PAML
    # with cwd set to it, so codeml's fixed-name outputs (2ML.dS, rst, rub)
    # cannot collide between workers. Note that concurrency multiplies peak
    # memory: every worker holds its own CDS index.
    total=0
    for aln in *.fas; do
        [ -e "\$aln" ] || continue
        total=\$((total + 1))
    done

    mkdir -p .ks_ok

    export KS_CDS="${cds_files}"
    export KS_ARGS="${args}"

    run_one() {
        aln="\$1"
        og="\${aln%%.*}"
        # \$KS_CDS and \$KS_ARGS are deliberately unquoted: both are
        # space-separated lists that must word-split into separate arguments.
        if ks.py \\
            --msa "\$aln" \\
            --cds \$KS_CDS \\
            --paralogs "\${og}.paralogs.txt" \\
            --tree "\${og}.tree" \\
            --orthogroup "\$og" \\
            --output "\${og}_ks.tsv" \\
            \$KS_ARGS; then
            if [ -s "\${og}_ks.tsv" ]; then
                # A marker file per success: counters set inside xargs workers
                # live in subshells and would be lost.
                touch ".ks_ok/\${og}"
            else
                echo "WARNING: ks.py completed but produced no output for \$og" >&2
            fi
        else
            echo "WARNING: ks.py failed for \$og" >&2
        fi
    }
    export -f run_one

    if [ "\${total}" -gt 0 ]; then
        printf '%s\\n' *.fas | xargs -r -I{} -P ${task.cpus} bash -c 'run_one "\$@"' _ {}
    fi

    successful=\$(find .ks_ok -type f | wc -l | tr -d ' ')
    rm -rf .ks_ok

    # A single orthogroup failing is expected and tolerated; only a batch in
    # which nothing succeeded is treated as a task failure.
    if [ "\${total}" -gt 0 ] && [ "\${successful}" -eq 0 ]; then
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
