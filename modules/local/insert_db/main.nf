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
        path(gbk_files)
  
  // global (single) files:
  path sketch_msh
  path db_in

  output:
  path "etd.db"       , emit: db
  path "register.log" , emit: log
  path "*_debug.txt", optional: true, emit: debug_files


  script:
  // build GBK argument for work directory files (for reading)
  def gbkArg = (gbk_files && gbk_files.size() > 0) ? 
    "--comp_gbk_files ${gbk_files.collect{ it.toString() }.join(' ')}" : ""

 // Build published path arguments (for DB storage)
  def amr_published = "${params.outdir}/amrfinder/${meta.id}.tsv"
  def contigs_published = "${params.outdir}/mobsuite/${meta.id}_results/contig_report.txt"
  def sketch_published = "${params.outdir}/mash/all_genomes.msh"
  def ice_published = ice_hits ? "${params.outdir}/filter/${meta.id}_ICEBERG_filtered_hits.tsv" : ""
  def phage_published = phage_coords ? "${params.outdir}/phispy/${meta.id}_phispy.tsv" : ""

  // Build published GBK paths
  def gbk_published = ""
  if (gbk_files && gbk_files.size() > 0) {
      def published_gbks = gbk_files.collect { "${params.outdir}/tncomp/${meta.id}_${it.name}" }
      gbk_published = "--comp_gbk_files_published ${published_gbks.join(' ')}"
    }

  """
  set -euo pipefail

  python3 ${projectDir}/bin/aggregate_output.py \\
    --db_path etd.db \\
    --fasta_name ${meta.id} \\
    --organism "${meta.organism}" \\
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
    ${gbkArg} \\
    ${gbk_published} \\

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
