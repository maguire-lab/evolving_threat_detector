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
  include { PARSE_GBK as PARSE_GBK_MAKE } from './modules/local/parsegbk/main'
  include { PARSE_GBK as PARSE_GBK_QUERY } from './modules/local/parsegbk/main'

  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      NAMED WORKFLOW FOR HELP MESSAGE
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  workflow ETD {
      if (params.help) {
          log.info """
          ==============================================
          ETD Pipeline - ETD (Evolving Threat Detector)
          ==============================================
          
          Usage:
            nextflow run main.nf --mode <MODE> [options]
          
          Modes:
            --mode make_db      Build database from reference genomes
            --mode query_db     Query existing database with new genomes
            --mode all          Run both workflows sequentially

          Required inputs:
            --input             Input samplesheet (used when input_make/query not specified)
            --input_make        Samplesheet for database building (make_db mode)
            --input_query       Samplesheet for database querying (query_db mode)

          Examples:
            # Build database only
            nextflow run main.nf --mode make_db --input reference_genomes.csv

            # Query existing database
            nextflow run main.nf --mode query_db --input query_genomes.csv

            # Run both sequentially
            nextflow run main.nf --mode all --input_make refs.csv --input_query queries.csv

            # Test the pipeline
            nextflow run main.nf -profile test

          For more information, visit: https://github.com/nf-core/etd
          ==============================================
          """.stripIndent()
          exit 0
      }

      main()
  }

  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      MAIN UNIFIED WORKFLOW
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  workflow {

      main:

      // Validate required parameters
      if (!params.mode) {
          error "Please specify --mode [make_db|query_db|all]."
      }

      if (!['make_db', 'query_db', 'all'].contains(params.mode)) {
          error "Invalid mode '${params.mode}'. Must be one of: make_db, query_db, all"
      }

      PIPELINE_INITIALISATION (
          params.version,
          params.validate_params,
          params.monochrome_logs,
          args,
          params.outdir,
          null
      )

      // Initialize variables to store outputs from MAKE_DB
      diamond_db_out = Channel.empty()
      reference_out = Channel.empty()
      etd_db_out = Channel.empty()

      //
      // WORKFLOW: Build database if requested
      //
      if (params.mode == 'make_db' || params.mode == 'all') {

          //log.info "Running MAKE_DB workflow..."

          // Determine input source for make_db
          def input_make_source = params.input_make ?: params.input
          if (!input_make_source) {
              error "No input specified for MAKE_DB. Use --input_make or --input"
          }

          // Validate input file exists
          if (!file(input_make_source).exists()) {
              error "Input file not found: ${input_make_source}"
          }

          //log.info "Using input for database building: ${input_make_source}"

          Channel
              .fromList(samplesheetToList(input_make_source,
  "${projectDir}/assets/schema_input.json"))
              .map { meta, gbk -> [meta, gbk] }
              .set { ch_samplesheet_make }

          ch_samplesheet_make | PARSE_GBK_MAKE

          PARSE_GBK_MAKE.out
              .map { meta, genome, protein, gff ->
                  tuple(meta, genome, protein, gff)
              }
              .set { ch_amrfinder_input_make }

          PARSE_GBK_MAKE.out
              .map { tuple ->
                  def (meta, genome, protein, gff) = tuple
                  [meta, genome]
              }
              .set { ch_genomes_make }

          PARSE_GBK_MAKE.out
              .map { tuple ->
                  def (meta, genome, protein, gff) = tuple
                  [meta, protein]
              }
              .set { ch_proteins_make }

          // Run the MAKE_DB subworkflow
          MAKE_DB_V1(ch_samplesheet_make, ch_amrfinder_input_make, ch_genomes_make, ch_proteins_make)

          // Store outputs for potential use by QUERY_DB
          diamond_db_out = MAKE_DB_V1.out.diamond_db ?: Channel.empty()
          reference_out = MAKE_DB_V1.out.sketch_reference ?: Channel.empty()
          etd_db_out = MAKE_DB_V1.out.updated_db ?: Channel.empty()

          //log.info "MAKE_DB workflow completed successfully"
      }

      //
      // WORKFLOW: Query database if requested
      //
      if (params.mode == 'query_db' || params.mode == 'all') {

          //log.info "Running QUERY_DB workflow..."

          // Determine input source for query_db
          def input_query_source = params.input_query ?: params.input
          if (!input_query_source) {
              error "No input specified for QUERY_DB. Use --input_query"
          }

          // Validate input file exists
          if (!file(input_query_source).exists()) {
              error "Input file not found: ${input_query_source}"
          }

          //log.info "Using input for database querying: ${input_query_source}"

          Channel
              .fromList(samplesheetToList(input_query_source,
  "${projectDir}/assets/schema_input.json"))
              .map { meta, gbk -> [meta, gbk] }
              .set { ch_samplesheet_query }

          ch_samplesheet_query | PARSE_GBK_QUERY

          PARSE_GBK_QUERY.out
              .map { meta, genome, protein, gff ->
                  tuple(meta, genome, protein, gff)
              }
              .set { ch_amrfinder_input_query }

          PARSE_GBK_QUERY.out
              .map { tuple ->
                  def (meta, genome, protein, gff) = tuple
                  [meta, genome]
              }
              .set { ch_genomes_query }

          PARSE_GBK_QUERY.out
              .map { tuple ->
                  def (meta, genome, protein, gff) = tuple
                  [meta, protein]
              }
              .set { ch_proteins_query }

          // Determine database sources based on mode
          if (params.mode == 'all') {
              //log.info "Using databases created by MAKE_DB workflow"
              diamond_db = diamond_db_out
              reference = reference_out
              etd_db_last = etd_db_out
          } else {
              // Use existing databases (query-only mode)
              //log.info "Loading existing databases from ${params.outdir}"

              diamond_db = Channel.fromPath("${params.outdir}/diamond/*.dmnd")
                  .collect()
                  .map { files ->
                      if (files.size() == 0) {
                          error "No .dmnd files found in ${params.outdir}/diamond/"
                      }
                      //log.info "Found diamond database: ${files[0]}"
                      return files[0]
                  }

              reference = Channel.value(file("${params.outdir}/mash/all_genomes.msh"))
                  .map { ref_file ->
                      if (!ref_file.exists()) {
                          error """Reference file not found: ${ref_file}. 
                                   Run MAKE_DB first or check --outdir path."""
                      }
                      //log.info "Found reference file: ${ref_file}"
                      return ref_file
                  }

              etd_db_last = Channel.value(file("${params.outdir}/insert/etd.db"))
                  .map { db_file ->
                      if (!db_file.exists()) {
                          error """Database file not found: ${db_file}. 
                                 Run MAKE_DB first or check --outdir path."""
                      }
                      //log.info "Found ETD database: ${db_file}"
                      return db_file
                  }
          }

          // Run the QUERY_DB subworkflow
          QUERY_DB_V1(ch_samplesheet_query, ch_amrfinder_input_query, ch_genomes_query,
  ch_proteins_query, reference, diamond_db, etd_db_last)

          //log.info "QUERY_DB workflow completed successfully"
      }

      //
      // WORKFLOW: Completion
      //
      //PIPELINE_COMPLETION (
          //params.email,
          //params.email_on_fail,
          //params.plaintext_email,
          //params.outdir,
          //params.monochrome_logs,
          //params.hook_url,
         // Channel.empty()
      //)

      //log.info " ETD Pipeline completed successfully!"
      //log.info "Results are available in: ${params.outdir}"
  }

  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      INDIVIDUAL ENTRY POINT WORKFLOWS (for development/testing)
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

  workflow MAKE_DB_ENTRY {

      main:

      //log.info "Running MAKE_DB_ENTRY workflow (development mode)"

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
          .map { meta, gbk -> [meta, gbk] }
          .set { ch_samplesheet }

      ch_samplesheet | PARSE_GBK

      PARSE_GBK.out
          .map { meta, genome, protein, gff ->
              tuple(meta, genome, protein, gff)
          }
          .set { ch_amrfinder_input }

      PARSE_GBK.out
          .map { tuple ->
              def (meta, genome, protein, gff) = tuple
              [meta, genome]
          }
          .set { ch_genomes }

      PARSE_GBK.out
          .map { tuple ->
              def (meta, genome, protein, gff) = tuple
              [meta, protein]
          }
          .set { ch_proteins }

      // Run the MAKE_DB subworkflow
      MAKE_DB_V1(ch_samplesheet, ch_amrfinder_input, ch_genomes, ch_proteins)
      

      emit:
      diamond_db = MAKE_DB_V1.out.diamond_db
      reference = MAKE_DB_V1.out.sketch_reference
      updated_db = MAKE_DB_V1.out.updated_db
  }

  workflow QUERY_DB_ENTRY {

      main:

      //log.info "Running QUERY_DB_ENTRY workflow (development mode)"

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
          .map { meta, gbk -> [meta, gbk] }
          .set { ch_samplesheet }

      ch_samplesheet | PARSE_GBK

      PARSE_GBK.out
          .map { meta, genome, protein, gff ->
              tuple(meta, genome, protein, gff)
          }
          .set { ch_amrfinder_input }

      PARSE_GBK.out
          .map { tuple ->
              def (meta, genome, protein, gff) = tuple
              [meta, genome]
          }
          .set { ch_genomes }

      PARSE_GBK.out
          .map { tuple ->
              def (meta, genome, protein, gff) = tuple
              [meta, protein]
          }
          .set { ch_proteins }

      // Load from published/cached results
      diamond_db = Channel.fromPath("${params.outdir}/diamond/*.dmnd")
          .collect()
          .map { files ->
              if (files.size() == 0) {
                  error "No .dmnd files found in ${params.outdir}/diamond/"
              }
              return files[0]
          }
      reference = Channel.value(file("${params.outdir}/mash/all_genomes.msh"))
      etd_db_last = Channel.value(file("${params.outdir}/insert/etd.db"))

      // Run the QUERY_DB subworkflow
      QUERY_DB_V1(ch_samplesheet, ch_amrfinder_input, ch_genomes, ch_proteins, reference,
  diamond_db, etd_db_last)

      emit:
      resistome_differences = QUERY_DB_V1.out.resistome_differences
  }

  // Backward compatibility - keep legacy names
  workflow MAKE_DB {
      MAKE_DB_ENTRY()
      
      emit:
      diamond_db = MAKE_DB_ENTRY.out.diamond_db
      reference = MAKE_DB_ENTRY.out.reference  
      updated_db = MAKE_DB_ENTRY.out.updated_db
  }

  workflow QUERY_DB {
      QUERY_DB_ENTRY()

      emit:
      resistome_differences = QUERY_DB_ENTRY.out.resistome_differences
  }

  workflow.onComplete {
      // Clean up etd.db from project root if it exists
      def root_db = file("etd.db")
      if (root_db.exists()) {
         root_db.delete()
         log.info "✓ Cleaned up etd.db from project root"
         log.info "✓ ETD Pipeline completed successfully!"
         log.info "✓ Results are available in: ${params.outdir}"
      
      }
}

  /*
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      THE END
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  */

