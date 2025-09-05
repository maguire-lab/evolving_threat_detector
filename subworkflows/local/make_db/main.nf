/*
 * Import required modules
 */
include { AMRFINDERPLUS_UPDATE      } from '../../../modules/nf-core/amrfinderplus/update'
include { AMRFINDERPLUS_RUN         } from '../../../modules/local/amrfinderplus/run'
include { MASH_SKETCH               } from '../../../modules/nf-core/mash/sketch'
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
include { INSERT_DB                 } from '../../../modules/local/insert_db'

workflow MAKE_DB {
    take:
    ch_samplesheet
    ch_amrfinder_input
    ch_genomes
    ch_proteins

    main:

    // Initialize database
    //etd_db_init = DB_INIT(Channel.value("etd.db"))

    // Insert genomes + organism into the DB
    //ch_db_insert = INSERT_GENOMES(ch_genomes, etd_db_init)

    // Update AMRFinderPlus database
    amrfinder_db = AMRFINDERPLUS_UPDATE()

    // Run MASH sketching
    MASH_SKETCH(ch_genomes)

    // Run MASH paste
    all_msh_list = MASH_SKETCH.out.mash.map { meta, msh -> msh }.collect()
    pasted = MASH_PASTE(all_msh_list)

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
    GET_ICEBERG()
    GET_ICEBERG.out.iceberg
        .set { ch_iceberg_db }

    // Step 2: Create DIAMOND database from ICEberg
    DIAMOND_MAKEDB(ch_iceberg_db)

    DIAMOND_MAKEDB.out.db.view { "DB: $it" }

    // Step 3: Run DIAMOND blastp with protein sequences from PARSE_GBK

    DIAMOND_BLASTP(
        ch_proteins,
        DIAMOND_MAKEDB.out.db,
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

    // Insert genomes into db
    // Step 1: Define short aliases for module outputs
    amr_ch    = AMRFINDERPLUS_RUN.out.report            
    mob_ch    = MOBSUITE_RECON.out.contig_report
        
    phage_ch  = PHISPY.out.coordinates  
    phage_ch.view { "phage_ch: $it" }
                
    ice_ch    = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
    ice_ch.view { "ice_ch: $it" }
  
    tncomp_all = TNCOMP_FINDER.out.gbk.groupTuple()
    tncomp_all.view { "After groupTuple: $it" }

    sketch_ref = MASH_PASTE.out.reference
    sketch_ref.view { "sketch_ref: $it" }

    db_initial = Channel.of(file(params.db_path ?: "etd.db")) 
    db_initial.view { "db_initial: $it" }  
 
    // Step 2: Join per genome by meta.id
    paired = amr_ch.combine(mob_ch, by: 0)              // [meta, amr_tsv, contigs_report]
    paired.view { "paired: $it" }

    paired2 = paired
       .combine(ice_ch, by: 0)                    // [meta, amr_tsv, contigs_report, ice_file_or_null]
    paired2.view { "paired2: $it" }
    paired3 = paired2
       .combine(phage_ch, by: 0)
    paired3.view { "paired3: $it" }

    paired4 = paired3                    // [meta, amr_tsv, contigs_report, ice_file_or_null, phage_file_or_null]
       .combine(tncomp_all, by: 0)                    // [meta, amr_tsv, contigs_report, ice_file_or_null, phage_file_or_null, gbk_files_list_or_null]

    paired4.view { "paired4: $it" }

    // Step 3: Shape the final per-genome bundle
    to_insert = paired4.map { meta, amr_tsv, contigs_report, ice_file, phage_file, gbk_files ->
    // Convert nulls to empty lists for optional inputs
          tuple(meta, amr_tsv, contigs_report,
           ice_file   ?: [], 
           phage_file ?: [], 
           gbk_files  ?: [])
}

    // Call INSERT_DB
    
     INSERT_DB(
        to_insert,
        sketch_ref,
        db_initial
    )

    emit:
    //db_etd               = etd_db_init.sqlite_db
    //updated_db           = ch_db_insert.sqlite_db
      mash_sketches        = MASH_SKETCH.out.mash
      sketch_reference     = MASH_PASTE.out.reference
      amr_reports          = AMRFINDERPLUS_RUN.out.report
      mobtyper_results     = MOBSUITE_RECON.out.contig_report
      integron_summaries   = INTEGRONFINDER.out.summary
      iceberg_hits         = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
      phispy_prophage_tsv  = PHISPY.out.prophage_tsv
      phispy_coordinates   = PHISPY.out.coordinates
      tncomp_gbk           = TNCOMP_FINDER.out.gbk
      updated_db           = INSERT_DB.out.db 
}
