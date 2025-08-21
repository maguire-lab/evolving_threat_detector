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

include { MAKE_DB                 } from './subworkflows/local/make_db'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_etd_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_etd_pipeline'
include { samplesheetToList       } from 'plugin/nf-schema'
include { PARSE_GBK               } from './modules/local/parsegbk/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DEVELOPMENT WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        null
    )

    Channel
        .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
        .map {
            meta, gbk ->
            [meta, gbk]
        }
        .set { ch_samplesheet }

    //ch_samplesheet.view { "$it" }

    ch_samplesheet | PARSE_GBK

    PARSE_GBK.out
        .map { 
            meta, genome, protein, gff ->
            tuple(meta, genome, protein, gff)
        }
        .set { ch_amrfinder_input }

    //ch_amrfinder_input.view { "$it" }

    PARSE_GBK.out
    .map { tuple ->
        def (meta, genome, protein, gff) = tuple
        [meta, genome]
    }
    .set { ch_genomes }

    //ch_genomes.view { "$it" }

    PARSE_GBK.out
    .map { tuple ->
        def (meta, genome, protein, gff) = tuple
        [meta, protein]
    }
    .set { ch_proteins }

    ch_proteins.view { "$it" }


    // Run the MAKE_DB subworkflow
    MAKE_DB(ch_samplesheet, ch_amrfinder_input, ch_genomes, ch_proteins)
    
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
