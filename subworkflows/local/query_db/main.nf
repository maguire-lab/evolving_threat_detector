/*
 * Import required modules
 */
include { AMRFINDERPLUS_UPDATE      } from '../../../modules/nf-core/amrfinderplus/update'
include { AMRFINDERPLUS_RUN         } from '../../../modules/local/amrfinderplus/run'
include { MASH_DIST                 } from '../../../modules/nf-core/mash/dist'
include { MASH_PASTE                } from '../../../modules/local/mash_paste'
include { MOBSUITE_RECON            } from '../../../modules/local/mobsuite/recon'
include { INTEGRONFINDER            } from '../../../modules/nf-core/integronfinder/main'
include { PARSE_GBK                 } from '../../../modules/local/parsegbk/main'
include { PHISPY                    } from '../../../modules/nf-core/phispy/main'
include { DB_INIT                   } from '../../../modules/local/db_init'
include { INSERT_GENOMES            } from '../../../modules/local/insert_genomes'
include { GET_ICEBERG               } from '../../../modules/local/iceberg/main'
include { DIAMOND_MAKEDB            } from '../../../modules/local/diamond/makedb'
include { DIAMOND_BLASTP            } from '../../../modules/local/diamond/blastp'
include { FILTER_DIAMOND_HITS       } from '../../../modules/local/filter_hits'
include { TNCOMP_FINDER             } from '../../../modules/local/tncomp_finder'
include { TN3_FINDER                } from '../../../modules/local/tn3finder'
include { QUERY_RESISTOME           } from '../../../modules/local/query_resistome'

workflow QUERY_DB {

    take:
    ch_samplesheet
    ch_amrfinder_input
    ch_genomes
    ch_proteins
    reference
    diamond_db
    etd_db_last

    main:

    // Initialize database
    //etd_db_init = DB_INIT(Channel.value("etd.db"))

    // Insert genomes + organism into the DB
    //ch_db_insert = INSERT_GENOMES(ch_genomes, etd_db_init)

    // Update AMRFinderPlus database
    amrfinder_db = AMRFINDERPLUS_UPDATE()

    // Run MASH distance estimation
    MASH_DIST(ch_genomes, reference)

    // Run MASH paste
    //all_msh_list = MASH_SKETCH.out.mash.map { meta, msh -> msh }.collect()
    //pasted = MASH_PASTE(all_msh_list)

    // Run AMRFinderPlus with updated DB
    AMRFINDERPLUS_RUN(ch_amrfinder_input, amrfinder_db[0])

    // Run MOB-suite
    MOBSUITE_RECON(ch_genomes)

    // Run IntegronFinder
    INTEGRONFINDER(ch_genomes)

    // Run Phispy
    PHISPY(ch_samplesheet)

   // Run ICEberg annotation

    // Step 1: Get the database
    //GET_ICEBERG()
    //GET_ICEBERG.out.iceberg
       // .set { ch_iceberg_db }

    // Step 2: Create DIAMOND database from ICEberg
    //DIAMOND_MAKEDB(ch_iceberg_db)

    //DIAMOND_MAKEDB.out.db.view { "DB: $it" }

    // Step 3: Run DIAMOND blastp with protein sequences from PARSE_GBK

    DIAMOND_BLASTP(
        ch_proteins,
        diamond_db,
        "txt",
        "qseqid sseqid pident slen qlen length mismatch gapopen qstart qend sstart send evalue bitscore"
    )
    DIAMOND_BLASTP.out.txt.view { "BlastP output channel: $it" }

    // Step 4: Filter the DIAMOND results
    ch_blast_results = DIAMOND_BLASTP.out.txt

    FILTER_DIAMOND_HITS(
        ch_blast_results,
        "ICEBERG",
        params.min_pident,
        params.min_alignment_length
    )

   // Run Tncompfinder
    TNCOMP_FINDER(ch_genomes)

   // Run Tn3finder
    TN3_FINDER(ch_genomes)

   // Query genomes in db
    // Step 1: Define short aliases for module outputs
    amr_ch    = AMRFINDERPLUS_RUN.out.report
    mob_ch    = MOBSUITE_RECON.out.contig_report

    phage_ch  = PHISPY.out.coordinates
    phage_ch.view { "phage_ch: $it" }

    ice_ch    = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
    ice_ch.view { "ice_ch: $it" }

    tncomp_all = TNCOMP_FINDER.out.gbk.groupTuple()
    tncomp_all.view { "After groupTuple: $it" }

    mash_dist = MASH_DIST.out.dist
    mash_dist.view { "mash_dist: $it" }

    // DEBUG: Count items in each channel
    amr_ch.count().view { "AMR reports count: $it" }
    mob_ch.count().view { "MOB reports count: $it" }
    phage_ch.count().view { "Phage coords count: $it" }
    ice_ch.count().view { "ICE hits count: $it" }
    tncomp_all.count().view { "TnComp grouped count: $it" }

    // DEBUG: View the meta IDs in each channel
    amr_ch.map { meta, files -> meta.id }.collect().view { "AMR IDs: $it" }
    mob_ch.map { meta, files -> meta.id }.collect().view { "MOB IDs: $it" }
    phage_ch.map { meta, files -> meta.id }.collect().view { "Phage IDs: $it" }
    ice_ch.map { meta, files -> meta.id }.collect().view { "ICE IDs: $it" }
    tncomp_all.map { meta, files -> meta.id }.collect().view { "TnComp IDs: $it" }


    etd_db_last =  Channel.fromPath("${params.outdir}/insert/etd.db")
    etd_db_last.view { "db_last: $it" }

    // Step 2: Join per genome by meta.id
    paired = amr_ch.join(mob_ch, by: 0, remainder: true)              // [meta, amr_tsv, contigs_report]
    paired.count().view { "After first join: $it genomes" }
    paired.view { "paired: $it" }

    paired2 = paired
       .join(ice_ch, by: 0, remainder: true)                    // [meta, amr_tsv, contigs_report, ice_file_or_null]
    paired2.count().view { "After second join: $it genomes" }
    paired2.view { "paired2: $it" }

    paired3 = paired2
       .join(phage_ch, by: 0, remainder: true)
    paired3.count().view { "After third join: $it genomes" }
    paired3.view { "paired3: $it" }

    paired4 = paired3                    // [meta, amr_tsv, contigs_report, ice_file_or_null, phage_file_or_null]
       .join(tncomp_all, by: 0, remainder: true)                    // [meta, amr_tsv, contigs_report, ice_file_or_null, phage_file_or_null, gbk_files_list_or_null]
    paired4.count().view { "After fourth join: $it genomes" }

    // Step 3: Shape the final per-genome bundle
    query_input = paired4
        .join(mash_dist, by: 0)
        .map { items ->
            def meta = items[0]
            def amr_tsv = items[1] ?: []
            def contigs_report = items[2] ?: []
            def ice_file = items[3] ?: []
            def phage_file = items[4] ?: []
            def gbk_files = items[5] ?: []
            def dist_file = items[6]
    
        tuple(meta, amr_tsv, contigs_report, ice_file, phage_file, gbk_files, dist_file)
}

    query_input.count().view { "Final query_input count: $it genomes" }
    query_input.view { "Query input structure: $it" }

     // Call QUERY_DB

     QUERY_RESISTOME(
        query_input,
        etd_db_last
    )


 emit:
    //db_etd               = etd_db_init.sqlite_db
    //updated_db           = ch_db_insert.sqlite_db
    //mash_sketches        = MASH_SKETCH.out.mash
    //sketch_reference     = MASH_PASTE.out.reference
      mash_distance        = MASH_DIST.out.dist
      amr_reports          = AMRFINDERPLUS_RUN.out.report
      mobtyper_results     = MOBSUITE_RECON.out.contig_report
      integron_summaries   = INTEGRONFINDER.out.summary
      iceberg_hits         = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
      phispy_prophage_tsv  = PHISPY.out.prophage_tsv
      phispy_coordinates   = PHISPY.out.coordinates
      tncomp_gbk           = TNCOMP_FINDER.out.gbk
}
