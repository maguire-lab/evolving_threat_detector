/*
 * Import required modules
 */
//include { AMRFINDERPLUS_UPDATE              } from '../../../modules/local/amrfinderplus/update'
include { AMRFINDERPLUS_RUN                 } from '../../../modules/local/amrfinderplus/run'
include { MASH_SKETCH                       } from '../../../modules/nf-core/mash/sketch'
include { MASH_PASTE_BATCH; MASH_PASTE_FINAL} from '../../../modules/local/mash_paste'
include { MOBSUITE_RECON                    } from '../../../modules/local/mobsuite/recon'
include { INTEGRONFINDER                    } from '../../../modules/nf-core/integronfinder/main'
include { PARSE_GBK                         } from '../../../modules/local/parsegbk/main'
include { PHISPY                            } from '../../../modules/nf-core/phispy/main'
include { GET_ICEBERG                       } from '../../../modules/local/iceberg/main'
include { DIAMOND_MAKEDB                    } from '../../../modules/local/diamond/makedb'
include { DIAMOND_BLASTP                    } from '../../../modules/local/diamond/blastp'
include { FILTER_DIAMOND_HITS               } from '../../../modules/local/filter_hits'
//include { TNCOMP_FINDER                     } from '../../../modules/local/tncomp_finder'
//include { TN3_FINDER                        } from '../../../modules/local/tn3finder'
include { INSERT_DB                         } from '../../../modules/local/insert_db'

workflow MAKE_DB {
    take:
    ch_samplesheet
    ch_amrfinder_input
    ch_genomes
    ch_proteins
    amrfinder_db

    main:

    // Update AMRFinderPlus database
    //amrfinder_db = AMRFINDERPLUS_UPDATE()

    // Run MASH sketching
    MASH_SKETCH(ch_genomes)

    // Run batched MASH paste
    // Collect all individual sketch files, then split into chunks of 50,000
    // This avoids the kernel max_map_count limit when mmap-ing 217K+ files
    // Works identically for small datasets — a single chunk is produced
    ch_all_msh = MASH_SKETCH.out.mash
        .map { meta, msh -> msh }
        .collect()

    ch_batches = ch_all_msh.flatMap { files ->
    // Normalise to a list (collect() may return a single item for 1 genome)
    def list = files instanceof List ? files : [files]
    def chunks = []
    // Split into sublists of up to 50,000 files each
    for (int i = 0; i < list.size(); i += 50000) {
        chunks << list.subList(i, Math.min(i + 50000, list.size()))
    }
    return chunks
}

    // Each chunk is pasted into an intermediate batch sketch
    MASH_PASTE_BATCH(ch_batches)

    // All batch sketches are collected and pasted into the final reference
    MASH_PASTE_FINAL(MASH_PASTE_BATCH.out.batch_msh.collect())

    // Run AMRFinderPlus with updated DB
    AMRFINDERPLUS_RUN(ch_amrfinder_input, amrfinder_db)

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

    // Step 3: Run DIAMOND blastp with protein sequences from PARSE_GBK

    DIAMOND_BLASTP(
        ch_proteins,
        DIAMOND_MAKEDB.out.db,
        "txt",
        "qseqid sseqid pident slen qlen length mismatch gapopen qstart qend sstart send evalue bitscore"
    )

    // Step 4: Filter the DIAMOND results
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

    // Insert genomes into db
    // Step 1: Define short aliases for module outputs
    amr_ch    = AMRFINDERPLUS_RUN.out.report            
    mob_ch    = MOBSUITE_RECON.out.contig_report
        
    phage_ch  = PHISPY.out.coordinates  
                
    ice_ch    = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
  
    //tncomp_all = TNCOMP_FINDER.out.report.groupTuple()

    //tn3_all = TN3_FINDER.out.report.groupTuple()

    tncomp_all = ch_genomes.map { meta, fasta -> tuple(meta, []) }

    tn3_all    = ch_genomes.map { meta, fasta -> tuple(meta, []) }

    sketch_ref = MASH_PASTE_FINAL.out.reference

    integron_ch = INTEGRONFINDER.out.integrons

    db_initial = Channel.of(file(params.db_path ?: "etd.db")) 
 
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

    // Join original GBK file from samplesheet (needed for ICE protein-contig mapping)
    paired7 = paired6
       .join(ch_samplesheet, by: 0, remainder: true)

    // Step 3: Shape the final per-genome bundle
    to_insert = paired7.map { items -> 
        def meta = items[0]
        def amr_tsv = items[1] ?: []
        def contigs_report = items[2] ?: []
        def ice_file = items[3] ?: []
        def phage_file = items[4] ?: []
        def txt_files_nested = items[5] ?: []
        def tn3_files_nested = items[6] ?: []
        def integron_file = items[7] ?: []
        def gbk_file = items[8] ?: []

     // Flatten the double-nested tncomp and tn3 gbk files
     def txt_files = txt_files_nested ? [txt_files_nested].flatten() : []
     
     def tn3_files = tn3_files_nested ? [tn3_files_nested].flatten() : []


       
        tuple(meta, amr_tsv, contigs_report, ice_file, phage_file, txt_files, tn3_files, integron_file, gbk_file)
}

    // Initialize empty database if it doesn't exist
    db_initial = file(params.db_path ?: "etd.db")
    if (!db_initial.exists()) {
        db_initial.text = ""
    }

    // Convert sketch_ref to a value channel so it can be reused
    sketch_value = sketch_ref.first()

    // Convert ICEberg raw FASTA to a value channel (global, shared across all genomes)
    iceberg_fasta_value = ch_iceberg_db.first()
    

    // Call INSERT_DB
    
     INSERT_DB(
        to_insert,
        sketch_value,
        db_initial,
        iceberg_fasta_value
    )

    emit:
      mash_sketches        = MASH_SKETCH.out.mash
      sketch_reference     = MASH_PASTE_FINAL.out.reference
      amr_reports          = AMRFINDERPLUS_RUN.out.report
      mobtyper_results     = MOBSUITE_RECON.out.contig_report
      integron_results     = INTEGRONFINDER.out.integrons
      diamond_db           = DIAMOND_MAKEDB.out.db
      iceberg_hits         = FILTER_DIAMOND_HITS.out.filtered_iceberg_hits
      phispy_prophage_tsv  = PHISPY.out.prophage_tsv
      phispy_coordinates   = PHISPY.out.coordinates
      //tncomp_report        = TNCOMP_FINDER.out.report
      //tn3_report           = TN3_FINDER.out.report
      updated_db           = INSERT_DB.out.db.last() 
      iceberg_fasta        = ch_iceberg_db
}	
