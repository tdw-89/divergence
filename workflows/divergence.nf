/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_divergence_pipeline'
include { ORTHOFINDER            } from '../modules/nf-core/orthofinder/main'
include { EXTRACT_PARALOGS       } from '../modules/local/extract_paralogs/main'
include { MAFFT_BATCH            } from '../modules/local/mafft_batch/main'
include { CODEML_BATCH           } from '../modules/local/codeml_batch/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Distribute items into batches of near-equal predicted cost.
//
// `.collate()` batches in channel order, and channel order here is filename
// order, and OrthoFinder numbers orthogroups by *descending* size. Collating
// therefore hands the first task the twenty largest orthogroups in the run --
// measured at 40% of a real dataset's total predicted cost in a single task,
// while the median task finished in under a minute and one ran past an eight
// hour wall clock. Dealing a cost-sorted list round-robin gives every batch a
// slice of the head and a slice of the tail, and leaves batch sizes within one
// item of each other. The cost model only has to get the *ordering* roughly
// right, so a rough proxy is enough.
def balancedBatches(List<List> costed, int size) {
    if (!costed) {
        return []
    }
    int bins = Math.max(1, (int) Math.ceil(costed.size() / (double) size))
    def batches = (0..<bins).collect { [] }
    costed.sort { a, b -> b[0] <=> a[0] }.eachWithIndex { entry, i ->
        batches[i % bins] << entry[1]
    }
    return batches
}

// Pairs the cross-checks will run on, straight from the paralog manifest:
// pairs are formed within each species, so this is the exact count, and it is
// what dominates a large orthogroup's cost.
def manifestPairCost(Path manifest) {
    def perSpecies = [:].withDefault { 0 }
    manifest.eachLine { line ->
        def parts = line.split('\t')
        if (parts.size() == 2 && parts[0] != 'species') {
            perSpecies[parts[0]] = perSpecies[parts[0]] + 1
        }
    }
    long pairs = 0
    perSpecies.each { _sp, n -> pairs += ((long) n) * (n - 1) / 2 }
    // Cross-checks stop at the cap, so beyond it extra pairs are free.
    if (params.max_pairwise_pairs && pairs > params.max_pairwise_pairs) {
        pairs = params.max_pairwise_pairs as long
    }
    // The M0 fit is the other half of the cost and grows with orthogroup size
    // rather than pair count. n^1.6 / 8 puts a 700-tip fit at roughly the same
    // number of cost units as the ~4 minutes it was measured to take.
    long genes = perSpecies.values().sum() ?: 0
    return pairs + (long) (Math.pow(genes, 1.6) / 8)
}


workflow DIVERGENCE {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    
    main:
    ch_versions = channel.empty()

    //
    // Collect all protein FASTAs from the samplesheet for OrthoFinder, and all
    // CDS FASTAs for the codon alignments further down.
    //
    ch_samplesheet
        .map { _meta, fasta, _cds -> fasta }
        .collect()
        .map { fastas -> [ [id: 'orthofinder'], fastas ] }
        .set { ch_fastas }

    ch_cds = ch_samplesheet
        .map { _meta, _fasta, cds -> cds }
        .collect()

    ch_prior_run = channel.of([ [:], [] ])

    //
    // MODULE: Run OrthoFinder
    //
    ORTHOFINDER(ch_fastas, ch_prior_run)

    //
    // MODULE: Select orthogroups containing target-species paralogs.
    // Emits the WHOLE orthogroup per selected group, plus the gene tree and a
    // manifest of the genes we actually want dS for.
    //
    EXTRACT_PARALOGS(ORTHOFINDER.out.orthofinder, params.target_species)

    //
    // CHANNEL PREPARATION: Batch orthogroup FASTAs for MAFFT
    //
    // Byte count stands in for cost here: MAFFT scales with the number of
    // sequences and their length, which is what the file size measures.
    ch_mafft_batches = EXTRACT_PARALOGS.out.fastas
        .flatten()
        .map { fasta -> [ fasta.size(), fasta ] }
        .toList()
        .flatMap { costed -> balancedBatches(costed, params.batch_size) }

    //
    // MODULE: Align each orthogroup (batched, L-INS-i)
    //
    MAFFT_BATCH(ch_mafft_batches)

    //
    // CHANNEL PREPARATION: Re-unite each alignment with its gene tree and
    // paralog manifest. All three share the <OG> filename stem; extraction
    // always writes a tree file (empty when OrthoFinder produced none) so this
    // join is total and no orthogroup is silently dropped.
    //
    def orthogroup_key = { file -> file.name.tokenize('.')[0] }

    ch_alignments = MAFFT_BATCH.out.fas.flatten().map { f -> [ orthogroup_key(f), f ] }
    ch_trees      = EXTRACT_PARALOGS.out.trees.flatten().map { f -> [ orthogroup_key(f), f ] }
    ch_manifests  = EXTRACT_PARALOGS.out.manifests.flatten().map { f -> [ orthogroup_key(f), f ] }

    ch_codeml_batches = ch_alignments
        .join(ch_trees)
        .join(ch_manifests)
        .map { _og, alignment, tree, manifest ->
            [ manifestPairCost(manifest), [ alignment, tree, manifest ] ]
        }
        .toList()
        .flatMap { costed -> balancedBatches(costed, params.batch_size) }
        .map { batch -> batch.flatten() }

    //
    // MODULE: Estimate pairwise dS between target-species paralogs (batched)
    //
    CODEML_BATCH(ch_codeml_batches, ch_cds)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }    

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'divergence_software_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    emit:
    orthofinder    = ORTHOFINDER.out.orthofinder     // channel: [ val(meta), path(orthofinder) ]
    summary        = EXTRACT_PARALOGS.out.summary    // channel: [ path(tsv) ]
    alignments     = MAFFT_BATCH.out.fas             // channel: [ path(fas) ]
    ks             = CODEML_BATCH.out.tsv            // channel: [ path(tsv) ]
    ds_trees       = CODEML_BATCH.out.ds_trees       // channel: [ path(nwk) ]
    dn_trees       = CODEML_BATCH.out.dn_trees       // channel: [ path(nwk) ]
    versions       = ch_versions                     // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
