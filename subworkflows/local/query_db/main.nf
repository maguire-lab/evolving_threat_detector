/*
 * Import required modules
 */
include { AMRFINDERPLUS_UPDATE      } from '../../../modules/nf-core/amrfinderplus/update'
include { AMRFINDERPLUS_RUN         } from '../../../modules/local/amrfinderplus/run'
include { MASH_DIST                 } from '../../../modules/nf-core/mash/dist'
include { MOBSUITE_RECON            } from '../../../modules/local/mobsuite/recon'
include { INTEGRONFINDER            } from '../../../modules/nf-core/integronfinder/main'
include { PARSE_GBK                 } from '../../../modules/local/parsegbk/main'
include { PHISPY                    } from '../../../modules/nf-core/phispy/main'
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

    // Update AMRFinderPlus database
    amrfinder_db = AMRFINDERPLUS_UPDATE()

    // Run MASH distance estimation
    MASH_DIST(ch_genomes, reference)

    // Run AMRFinderPlus with updated DB
    AMRFINDERPLUS_RUN(ch_amrfinder_input, amrfinder_db[0])

    // Run MOB-suite
    MOBSUITE_RECON(ch_genomes)

    // Run IntegronFinder
    INTEGRONFINDER(ch_genomes)

    // Run Phispy
    PHISPY(ch_samplesheet)

    // Run ICEberg annotation

    // Step 1: Run DIAMOND blastp with existing diamond ice db

    DIAMOND_BLASTP(
        ch_proteins,
        diamond_db,
        "txt",
        "qseqid sseqid pident slen qlen length mismatch gapopen qstart qend sstart send evalue bitscore"
    )
    //DIAMOND_BLASTP.out.txt.view { "BlastP output channel: $it" }

    // Step 2: Filter the DIAMOND results
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
    //phage_ch.view { "phage_ch: $it" }

    ice_ch    = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
    //ice_ch.view { "ice_ch: $it" }

    tncomp_all = TNCOMP_FINDER.out.gbk.groupTuple()
    //tncomp_all.view { "After groupTuple: $it" }

    tn3_all = TN3_FINDER.out.gbk.groupTuple()

    mash_dist = MASH_DIST.out.dist
    //mash_dist.view { "mash_dist: $it" }

    // DEBUG: Count items in each channel
    //amr_ch.count().view { "AMR reports count: $it" }
    //mob_ch.count().view { "MOB reports count: $it" }
    //phage_ch.count().view { "Phage coords count: $it" }
    //ice_ch.count().view { "ICE hits count: $it" }
    //tncomp_all.count().view { "TnComp grouped count: $it" }

    // DEBUG: View the meta IDs in each channel
    //amr_ch.map { meta, files -> meta.id }.collect().view { "AMR IDs: $it" }
    //mob_ch.map { meta, files -> meta.id }.collect().view { "MOB IDs: $it" }
    //phage_ch.map { meta, files -> meta.id }.collect().view { "Phage IDs: $it" }
    //ice_ch.map { meta, files -> meta.id }.collect().view { "ICE IDs: $it" }
    //tncomp_all.map { meta, files -> meta.id }.collect().view { "TnComp IDs: $it" }


    //etd_db_last =  Channel.value(file("${params.outdir}/insert/etd.db"))
    //etd_db_last.view { "db_last: $it" }

    // Step 2: Join per genome by meta.id
    paired = amr_ch.join(mob_ch, by: 0, remainder: true)
    //paired.count().view { "After first join: $it genomes" }
    //paired.view { "paired: $it" }

    paired2 = paired
       .join(ice_ch, by: 0, remainder: true)
    //paired2.count().view { "After second join: $it genomes" }
    //paired2.view { "paired2: $it" }

    paired3 = paired2
       .join(phage_ch, by: 0, remainder: true)
    //paired3.count().view { "After third join: $it genomes" }
    //paired3.view { "paired3: $it" }

    paired4 = paired3
       .join(tncomp_all, by: 0, remainder: true)
    //paired4.count().view { "After fourth join: $it genomes" }

    paired5 = paired4
       .join(tn3_all, by: 0, remainder: true)

    // Step 3: Shape the final per-genome bundle
    query_input = paired5
        .join(mash_dist, by: 0, remainder: true)
        .map { items ->
            def meta = items[0]
            def amr_tsv = items[1] ?: []
            def contigs_report = items[2] ?: []
            def ice_file = items[3] ?: []
            def phage_file = items[4] ?: []
            def gbk_files_nested = items[5] ?: []
            def tn3_files_nested = items[6] ?: []
            def dist_file = items[7]
      
     def gbk_files = gbk_files_nested ? [gbk_files_nested].flatten() : []
     def tn3_files = tn3_files_nested ? [tn3_files_nested].flatten() : []

        tuple(meta, amr_tsv, contigs_report, ice_file, phage_file, gbk_files, tn3_files, dist_file)
}

    //query_input.count().view { "Final query_input count: $it genomes" }
    //query_input.view { "Query input structure: $it" }

    // Call QUERY_DB

    QUERY_RESISTOME(
       query_input,
       etd_db_last
    )


 emit:
      mash_distance        = MASH_DIST.out.dist
      amr_reports          = AMRFINDERPLUS_RUN.out.report
      mobtyper_results     = MOBSUITE_RECON.out.contig_report
      integron_summaries   = INTEGRONFINDER.out.summary
      iceberg_hits         = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
      phispy_prophage_tsv  = PHISPY.out.prophage_tsv
      phispy_coordinates   = PHISPY.out.coordinates
      tncomp_gbk           = TNCOMP_FINDER.out.gbk
      tn3_gbk 		   = TN3_FINDER.out.gbk
}
