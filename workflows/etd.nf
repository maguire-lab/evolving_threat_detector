/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { MAKE_DB } from '../subworkflows/local/make_db'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ETD {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:

    // optional: track versions of my tools further down the line
    ch_versions = Channel.empty()

    // Collate and save software versions

    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_etd_software_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    // Prepare the (meta, fasta) tuples expected by modules
    ch_genomes = ch_samplesheet.map { file -> tuple([id: file.basename], file) }

    // Run the make_db workflow
    MAKE_DB(ch_genomes)


    emit:
    mash_sketches	= MAKE_DB.out.mash_sketches
    amr_reports 	= MAKE_DB.out.amr_reports
    mob_reports 	= MAKE_DB.out.mob_recon_results
    integron_results = MAKE_DB.out.integron_summaries
    versions       	= ch_versions

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
