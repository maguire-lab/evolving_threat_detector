/*
 * Import required modules
 */
include { AMRFINDERPLUS_UPDATE } from '../../../modules/nf-core/amrfinderplus/update'
include { AMRFINDERPLUS_RUN    } from '../../../modules/nf-core/amrfinderplus/run'
include { MASH_SKETCH          } from '../../../modules/nf-core/mash/sketch'
include { MOBSUITE_RECON       } from '../../../modules/nf-core/mobsuite/recon'
include { INTEGRONFINDER       } from '../../../modules/nf-core/integronfinder/main'
include { PROKKA               } from '../../../modules/nf-core/prokka/main'
include { PHISPY               } from '../../../modules/nf-core/phispy/main'
include { DB_INIT              } from '../../../modules/local/db_init'

workflow MAKE_DB {
    take:
    ch_genomes

    main:

    // Initialize database
    etd_db = DB_INIT(Channel.value("etd.db"))

    // Update AMRFinderPlus database
    amrfinder_db = AMRFINDERPLUS_UPDATE()

    // Run MASH sketching
    MASH_SKETCH(ch_genomes)

    // Run AMRFinderPlus with updated DB
    AMRFINDERPLUS_RUN(ch_genomes, amrfinder_db[0])

    // Run MOB-suite
    MOBSUITE_RECON(ch_genomes)

    // Run IntegronFinder
    INTEGRONFINDER(ch_genomes)

    // Generate .gbk file for phispy
    prokka_out = PROKKA(ch_genomes)
    
     //ch_gbk = prokka_out.gbk
         //.map { genomeID, gbk -> tuple(genomeID, gbk) }
    

    // Run Phispy
    //PHISPY(ch_gbk).ext { genomeID, gbk -> [ prefix: "${genomeID}_phispy" ] }

    emit:
    db_etd               = etd_db.sqlite_db
    //mash_sketches        = MASH_SKETCH.out.mash
    //amr_reports        = AMRFINDERPLUS_RUN.out.report
    //mobtyper_results   = MOBSUITE_RECON.out.contig_report
    //integron_summaries = INTEGRONFINDER.out.summary
    //ch_gbk_wprefix       = PROKKA.out.gbk
}
