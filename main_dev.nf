#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/etd
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/etd
    Website: https://nf-co.re/etd
    Slack  : https://nfcore.slack.com/channels/etd
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MAKE_DB } from './subworkflows/local/make_db'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_etd_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_etd_pipeline'
include { samplesheetToList       } from 'plugin/nf-schema'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DEVELOPMENT WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    // Create a channel from user-supplied input FASTA files
    //ch_input = Channel.from(params.input)
        //.map { row -> tuple(id: row.sample, file(row.fasta)) }

    //
    // SUBWORKFLOW: Run initialisation tasks
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        null
    )

    //PIPELINE_INITIALISATION.out.samplesheet.view { "SAMPLE: $it" }

   // PIPELINE_INITIALISATION.out.samplesheet
       // .map { meta, fastas -> tuple([id:meta.id], fastas[0]) }
      // .set { ch_input }

    //ch_input = Channel
        //.fromPath(params.input, checkIfExists: true)
        //.map { file -> tuple(file.baseName], file) }
        //.set { ch_input }

   // Channel
   //     .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
   //     .view()
   //     .map {
   //         meta, fasta ->
   //          return [ meta.id, meta + [ single_end:true ], [ fasta ] ]
   //     }
   //     .groupTuple()
   //     .map {
   //         meta, fasta ->
   //             return [ meta, fasta.flatten() ]
   //     }
   //     .set { ch_samplesheet }
    Channel
        .fromPath(params.input, checkIfExists: true)
        .map { file -> tuple(file.baseName, file) }
        .set { genome }

     genome.view { "$it" }

    // Run the MAKE_DB subworkflow with real data
    MAKE_DB(genome)
    
    // SUBWORKFLOW: Run completion tasks
   // PIPELINE_COMPLETION (
   //     params.email,
   //     params.email_on_fail,
   //     params.plaintext_email,
   //     params.outdir,
   //    params.monochrome_logs,
   //    params.hook_url,
   //     Channel.empty()
    //)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
