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
include { MAFFT_ALIGN            } from '../modules/nf-core/mafft/align/main'

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
    // Collect all FASTA files from the samplesheet for OrthoFinder
    //
    ch_samplesheet
        .map { meta, fasta -> fasta }
        .collect()
        .map { fastas -> [ [id: 'orthofinder'], fastas ] }
        .set { ch_fastas }

    ch_prior_run = channel.of([ [:], [] ])

    //
    // MODULE: Run OrthoFinder
    //
    ORTHOFINDER(ch_fastas, ch_prior_run)

    //
    // MODULE: Extract Paralogs
    // Pass the OrthoFinder output directory and your target species parameter
    //
    EXTRACT_PARALOGS(ORTHOFINDER.out.orthofinder, params.target_species)

    //
    // CHANNEL PREPARATION: Setup for MAFFT
    // Take the bulk output of EXTRACT_PARALOGS, flatten it into individual files,
    // and attach a unique 'meta' map to each so nf-core modules can process them.
    //
    ch_for_mafft = EXTRACT_PARALOGS.out.unaligned_fastas
        .flatten()
        .map { file -> 
            def meta = [id: file.baseName] 
            return [meta, file]
        }

    //
    // MODULE: Align paralog sequences
    //
    ch_dummy = channel.value([ [:], [] ])
    MAFFT_ALIGN(
        ch_for_mafft,
        ch_dummy,
        ch_dummy,
        ch_dummy,
        ch_dummy,
        ch_dummy,
        ch_dummy
    )

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
    alignments     = MAFFT_ALIGN.out.fas             // channel: [ val(meta), path(fas) ]
    versions       = ch_versions                     // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
