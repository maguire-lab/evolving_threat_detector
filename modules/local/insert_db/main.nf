process INSERT_DB {
  tag "${meta.id}"
  label 'process_low'

  container "quay.io/biocontainers/biopython:1.84"

 // publishDir "${params.outdir}/database", mode: 'copy'

  maxForks 1

  input:
  // per-genome bundle prepared by joins/mapping
  tuple val(meta), path(amr_tsv), path(contigs_report),
        path(ice_hits),
        path(phage_coords),
        path(txt_files),
        path(tn3_files),
        path(integron_file),
        path(gbk_file)
  
  // global (single) files:
  path sketch_msh
  path db_in
  path iceberg_fasta
  path geometry_tsvs

  output:
  path "etd.db"       , emit: db
  path "register.log" , emit: log
  path "*_debug.txt", optional: true, emit: debug_files


  script:
  // Handle organism
  def organism_value = meta.organism instanceof List ? "" : (meta.organism ?: "")

  // build argument for tncomp work directory files (for reading)
  def compArg = (txt_files && txt_files.size() > 0) ? 
    "--comp_txt_files ${txt_files.collect{ it.toString() }.join(' ')}" : ""

 // build TN3 argument for tn3 work directory files
  def tn3Arg = (tn3_files && tn3_files.size() > 0) ?
    "--tn3_txt_files ${tn3_files.collect{ it.toString() }.join(' ')}" : ""

  def integronArg = (integron_file && integron_file.size() > 0) ?
    "--integron_file ${integron_file}" : ""

 // Build published path arguments (for DB storage)
  def amr_published = "${params.outdir}/amrfinder/${meta.id}.tsv"
  def contigs_published = "${params.outdir}/mobsuite/${meta.id}_results/contig_report.txt"
  def sketch_published = "${params.outdir}/mash/all_genomes.msh"
  def ice_published = ice_hits ? "${params.outdir}/filter/${meta.id}_ICEBERG_filtered_hits.tsv" : ""
  def phage_published = phage_coords ? "${params.outdir}/phispy/${meta.id}_phispy.tsv" : ""

  // Build published transposons and integron paths

  def comp_published = ""
  if (txt_files && txt_files.size() > 0) {
      def published_txts = txt_files.collect { "${params.outdir}/tncomp/${meta.id}_${it.name}" }
      comp_published = "--comp_txt_files_published ${published_txts.join(' ')}"
  }

  def tn3_published = ""
  if (tn3_files && tn3_files.size() > 0) {
      def published_tn3s = tn3_files.collect { "${params.outdir}/tn3/${meta.id}_${it.name}" }
      tn3_published = "--tn3_txt_files_published ${published_tn3s.join(' ')}"
  }

  def integron_published = (integron_file && integron_file.size() > 0) ?
    "${params.outdir}/integronfinder/${meta.id}.integrons" : ""


  """
  set -euo pipefail

  # Seed from provided DB (previous task's output or initial DB)
  #if [ -s "seed.db" ]; then
      #cp "seed.db" etd.db
  #fi

  python3 ${projectDir}/bin/aggregate_output.py \\
    --db_path etd.db \\
    --fasta_name ${meta.id} \\
    ${organism_value ? "--organism \"${organism_value}\"" : ""} \\
    --sketch_path "${sketch_msh}" \\
    --sketch_path_published "${sketch_published}" \\
    --amrfinder_output "${amr_tsv}" \\
    --amrfinder_output_published "${amr_published}" \\
    --contigs_report_path "${contigs_report}" \\
    --contigs_report_path_published "${contigs_published}" \\
    ${ice_hits && ice_hits.size() > 0 ? "--filtered_hits_report_path ${ice_hits}" : ""} \\
    ${ice_hits && ice_hits.size() > 0 ? "--filtered_hits_report_path_published ${ice_published}" : ""} \\
    ${phage_coords && phage_coords.size() > 0 ? "--phage_report_path ${phage_coords}" : ""} \\
    ${phage_coords && phage_coords.size() > 0 ? "--phage_report_path_published ${phage_published}" : ""} \\
    ${compArg} \\
    ${comp_published} \\
    ${tn3Arg} \\
    ${tn3_published}\\
    ${integronArg} \\
    ${integron_published ? "--integron_file_published ${integron_published}" : ""} \\
    ${gbk_file && gbk_file.size() > 0 ? "--gbk_path ${gbk_file}" : ""} \\
    --iceberg_fasta ${iceberg_fasta} \\
    --max_distance ${params.max_distance} \\
    2>&1 | tee register.log

  # Capture any debug files created
  ls *_debug.txt 2>/dev/null || true
  """

  stub:
  """
  touch etd.db
  touch register.log
  touch ${meta.id}_debug.txt
  """
}

process INSERT_DB_BATCH {
  label 'process_medium'

  container "quay.io/biocontainers/biopython:1.84"

  //publishDir "${params.outdir}/insert", mode: 'copy'

  input:
  path manifest
  path sketch_msh
  path db_in, stageAs: 'seed.db'
  path iceberg_fasta
  path geometry_tsvs

  output:
  path "etd.db",       emit: db
  path "register.log", emit: log

  script:
  def iceberg_arg = iceberg_fasta.name != 'NO_ICEBERG' ? "--iceberg_fasta ${iceberg_fasta}" : ""
  """
  set -euo pipefail

  # Start from seed database if it has content; otherwise start fresh
  if [ -s seed.db ]; then
    cp seed.db etd.db
  fi

  python3 ${projectDir}/bin/batch_insert.py \\
    --manifest ${manifest} \\
    --db_path etd.db \\
    --sketch_path_published "${params.outdir}/mash/all_genomes.msh" \\
    ${iceberg_arg} \\
    --max_distance ${params.max_distance} \\
    --outdir ${params.outdir} \\
    --commit_every 500 \\
    2>&1 | tee register.log
  """

  stub:
  """
  touch etd.db
  touch register.log
  """
}
