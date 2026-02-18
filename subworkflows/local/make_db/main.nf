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

    //DIAMOND_MAKEDB.out.db.view { "DB: $it" }

    // Step 3: Run DIAMOND blastp with protein sequences from PARSE_GBK

    DIAMOND_BLASTP(
        ch_proteins,
        DIAMOND_MAKEDB.out.db,
        "txt",
        "qseqid sseqid pident slen qlen length mismatch gapopen qstart qend sstart send evalue bitscore"
    )
    //DIAMOND_BLASTP.out.txt.view { "BlastP output channel: $it" }

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
    //phage_ch.view { "phage_ch: $it" }
                
    ice_ch    = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
    //ice_ch.view { "ice_ch: $it" }
  
    tncomp_all = TNCOMP_FINDER.out.gbk.groupTuple()
    //tncomp_all.view { "After groupTuple: $it" }

    tn3_all = TN3_FINDER.out.gbk.groupTuple()

    sketch_ref = MASH_PASTE.out.reference
    //sketch_ref.view { "sketch_ref: $it" }

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


    db_initial = Channel.of(file(params.db_path ?: "etd.db")) 
    //db_initial.view { "db_initial: $it" }  
 
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
    //paired4.view { "paired4: $it" }

    paired5 = paired4
       .join(tn3_all, by: 0, remainder: true)

    // Step 3: Shape the final per-genome bundle
    to_insert = paired5.map { items -> 
        def meta = items[0]
        def amr_tsv = items[1] ?: []
        def contigs_report = items[2] ?: []
        def ice_file = items[3] ?: []
        def phage_file = items[4] ?: []
        def gbk_files_nested = items[5] ?: []
        def tn3_files_nested = items[6] ?: []

     // Flatten the double-nested tncomp and tn3 gbk files
     def gbk_files = gbk_files_nested ? [gbk_files_nested].flatten() : []
     
     def tn3_files = tn3_files_nested ? [tn3_files_nested].flatten() : []


       
        tuple(meta, amr_tsv, contigs_report, ice_file, phage_file, gbk_files, tn3_files)
}

    //to_insert.count().view { "Final to_insert count: $it genomes" }
    //to_insert.view {"Final to insert: $it" }

    // Initialize empty database if it doesn't exist
    db_initial = file(params.db_path ?: "etd.db")
    if (!db_initial.exists()) {
        db_initial.text = ""
    }
    //db_initial.view { "db_initial: $it" }

    // Convert sketch_ref to a value channel so it can be reused
    sketch_value = sketch_ref.first()
    

    // Call INSERT_DB
    
     INSERT_DB(
        to_insert,
        sketch_value,
        db_initial
    )

    emit:
      mash_sketches        = MASH_SKETCH.out.mash
      sketch_reference     = MASH_PASTE.out.reference
      amr_reports          = AMRFINDERPLUS_RUN.out.report
      mobtyper_results     = MOBSUITE_RECON.out.contig_report
      integron_summaries   = INTEGRONFINDER.out.summary
      diamond_db           = DIAMOND_MAKEDB.out.db
      iceberg_hits         = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
      phispy_prophage_tsv  = PHISPY.out.prophage_tsv
      phispy_coordinates   = PHISPY.out.coordinates
      tncomp_gbk           = TNCOMP_FINDER.out.gbk
      tn3_gbk		   = TN3_FINDER.out.gbk
      updated_db           = INSERT_DB.out.db.last() 
}	
