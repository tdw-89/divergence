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
    ch_mafft_batches = EXTRACT_PARALOGS.out.fastas
        .flatten()
        .collate(params.batch_size)

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
        .map { _og, alignment, tree, manifest -> [ alignment, tree, manifest ] }
        .collate(params.batch_size)
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
