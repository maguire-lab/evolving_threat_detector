/*
 * Import required modules
 */
//include { AMRFINDERPLUS_UPDATE      } from '../../../modules/local/amrfinderplus/update'
include { AMRFINDERPLUS_RUN         } from '../../../modules/local/amrfinderplus/run'
include { MASH_DIST                 } from '../../../modules/nf-core/mash/dist'
include { MOBSUITE_RECON            } from '../../../modules/local/mobsuite/recon'
include { INTEGRONFINDER            } from '../../../modules/nf-core/integronfinder/main'
include { PARSE_GBK                 } from '../../../modules/local/parsegbk/main'
include { PHISPY                    } from '../../../modules/nf-core/phispy/main'
include { DIAMOND_BLASTP            } from '../../../modules/local/diamond/blastp'
include { FILTER_DIAMOND_HITS       } from '../../../modules/local/filter_hits'
//include { TNCOMP_FINDER             } from '../../../modules/local/tncomp_finder'
//include { TN3_FINDER                } from '../../../modules/local/tn3finder'
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
    iceberg_fasta
    amrfinder_db

    main:

    // Update AMRFinderPlus database
    //amrfinder_db = AMRFINDERPLUS_UPDATE()

    // Run MASH distance estimation
    MASH_DIST(ch_genomes, reference)

    // Run AMRFinderPlus with updated DB
    AMRFINDERPLUS_RUN(ch_amrfinder_input, amrfinder_db)

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
    //TNCOMP_FINDER(ch_genomes)

   // Run Tn3finder
    //TN3_FINDER(ch_genomes)

    // Query genomes in db
    // Step 1: Define short aliases for module outputs
    amr_ch    = AMRFINDERPLUS_RUN.out.report
    mob_ch    = MOBSUITE_RECON.out.contig_report

    phage_ch  = PHISPY.out.coordinates

    ice_ch    = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits

    //tncomp_all = TNCOMP_FINDER.out.report.groupTuple()

    //tn3_all = TN3_FINDER.out.report.groupTuple()

    tncomp_all = ch_genomes.map { meta, fasta -> tuple(meta, []) }

    tn3_all = ch_genomes.map { meta, fasta -> tuple(meta, []) }

    mash_dist = MASH_DIST.out.dist

    integron_ch = INTEGRONFINDER.out.integrons


    //etd_db_last =  Channel.value(file("${params.outdir}/insert/etd.db"))
    //etd_db_last.view { "db_last: $it" }

    // Step 2: Join per genome by meta.id
    paired = amr_ch.join(mob_ch, by: 0, remainder: true)

    paired2 = paired
       .join(ice_ch, by: 0, remainder: true)

    paired3 = paired2
       .join(phage_ch, by: 0, remainder: true)

    paired4 = paired3
       .join(tncomp_all, by: 0, remainder: true)

    paired5 = paired4
       .join(tn3_all, by: 0, remainder: true)

    paired6 = paired5
       .join(integron_ch, by: 0, remainder: true)

    paired7 = paired6
       .join(ch_samplesheet, by: 0, remainder: true)


    // Step 3: Shape the final per-genome bundle
    query_input = paired7
        .join(mash_dist, by: 0, remainder: true)
        .map { items ->
            def meta = items[0]
            def amr_tsv = items[1] ?: []
            def contigs_report = items[2] ?: []
            def ice_file = items[3] ?: []
            def phage_file = items[4] ?: []
            def txt_files_nested = items[5] ?: []
            def tn3_files_nested = items[6] ?: []
            def integron_file = items[7] ?: []
            def gbk_file      = items[8] ?: []
            def dist_file = items[9]
      
     def txt_files = txt_files_nested ? [txt_files_nested].flatten() : []
     def tn3_files = tn3_files_nested ? [tn3_files_nested].flatten() : []

        tuple(meta, amr_tsv, contigs_report, ice_file, phage_file, txt_files, tn3_files, integron_file, dist_file, gbk_file)
}

    // Call QUERY_DB

    QUERY_RESISTOME(
       query_input,
       etd_db_last,
       iceberg_fasta
    )


 emit:
      mash_distance        = MASH_DIST.out.dist
      amr_reports          = AMRFINDERPLUS_RUN.out.report
      mobtyper_results     = MOBSUITE_RECON.out.contig_report
      integron_rsults      = INTEGRONFINDER.out.integrons
      iceberg_hits         = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
      phispy_prophage_tsv  = PHISPY.out.prophage_tsv
      phispy_coordinates   = PHISPY.out.coordinates
      //tncomp_report        = TNCOMP_FINDER.out.report
      //tn3_report           = TN3_FINDER.out.report
}
