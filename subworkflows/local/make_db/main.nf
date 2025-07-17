/*
 * Import required modules
 */
include { AMRFINDERPLUS_UPDATE } from '../../../modules/nf-core/amrfinderplus/update'
include { AMRFINDERPLUS_RUN    } from '../../../modules/local/amrfinderplus/run'
include { MASH_SKETCH          } from '../../../modules/nf-core/mash/sketch'
include { MOBSUITE_RECON       } from '../../../modules/local/mobsuite/recon'
include { INTEGRONFINDER       } from '../../../modules/nf-core/integronfinder/main'
include { PARSE_GBK            } from '../../../modules/local/parsegbk/main'
include { PHISPY               } from '../../../modules/nf-core/phispy/main'
include { DB_INIT              } from '../../../modules/local/db_init'
include { INSERT_GENOMES       } from '../../../modules/local/insert_genomes'

workflow MAKE_DB {
    take:
    ch_samplesheet
    ch_amrfinder_input
    ch_genomes

    main:

    // Initialize database
    //etd_db_init = DB_INIT(Channel.value("etd.db"))

    // Insert genomes + organism into the DB
    //ch_db_insert = INSERT_GENOMES(ch_genomes, etd_db_init)

    // Update AMRFinderPlus database
    amrfinder_db = AMRFINDERPLUS_UPDATE()

    // Run MASH sketching
    MASH_SKETCH(ch_genomes)

    // Run AMRFinderPlus with updated DB
    AMRFINDERPLUS_RUN(ch_amrfinder_input, amrfinder_db[0])

    // Run MOB-suite
    MOBSUITE_RECON(ch_genomes)

    // Run IntegronFinder
    INTEGRONFINDER(ch_genomes)
    
    // Run Phispy
    PHISPY(ch_samplesheet)

    emit:
    //db_etd               = etd_db_init.sqlite_db
    //updated_db           = ch_db_insert.sqlite_db
    mash_sketches        = MASH_SKETCH.out.mash
    amr_reports        = AMRFINDERPLUS_RUN.out.report
    mobtyper_results   = MOBSUITE_RECON.out.contig_report
    integron_summaries = INTEGRONFINDER.out.summary
}
