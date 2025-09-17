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

include { MAKE_DB as MAKE_DB_V1   } from './subworkflows/local/make_db'
include { QUERY_DB as QUERY_DB_V1 } from './subworkflows/local/query_db'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_etd_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_etd_pipeline'
include { samplesheetToList       } from 'plugin/nf-schema'
include { PARSE_GBK               } from './modules/local/parsegbk/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN DEVELOPMENT WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MAKE_DB {

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
    MAKE_DB_V1(ch_samplesheet, ch_amrfinder_input, ch_genomes, ch_proteins)
    
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

workflow QUERY_DB {

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

    // Load from published/cached results
    diamond_db = Channel.fromPath("${params.outdir}/diamond/*.dmnd")
    reference = Channel.fromPath("${params.outdir}/mash/all_genomes.msh")
    etd_db_last =  Channel.fromPath("${params.outdir}/insert/etd.db")

    // Run the QUERY_DB subworkflow
    QUERY_DB_V1(ch_samplesheet, ch_amrfinder_input, ch_genomes, ch_proteins, reference, diamond_db, etd_db_last)

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
