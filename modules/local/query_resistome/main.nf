process QUERY_RESISTOME {
  tag "${meta.id}"
  label 'process_medium'

  container "quay.io/precious/biopython-polars@sha256:056f5ae506d293c3bbc8de856e8b1cb5eb74c6179a8aa660ceb85b6828393bc9"

  publishDir "${params.outdir}/resistome_analysis", mode: 'copy'

  input:
  // per-genome bundle prepared by joins/mapping
  tuple val(meta), path(amr_tsv), path(contigs_report),
        path(ice_hits),
        path(phage_coords),
        path(txt_files),
        path(tn3_files),
        path(integron_file),
        path(mash_dist_output)
  
  // global (single) files:
  path db_in

  output:
  path "${meta.id}_resistome_differences.*"       , emit: differences
  path "*.log"                , emit: logs, optional: true
  
  when:
  task.ext.when == null || task.ext.when


  script:
  // Handle organism
  def organism_value = meta.organism instanceof List ? "" : (meta.organism ?: "")

  // build transposons and integron finders argument for directory files

  def compArg = (txt_files && txt_files.size() > 0) ? 
    "--comp_txt_files ${txt_files.collect{ it.toString() }.join(' ')}" : ""

  def tn3Arg = (tn3_files && tn3_files.size() > 0) ?
    "--tn3_txt_files ${tn3_files.collect{ it.toString() }.join(' ')}" : ""

  def integronArg = (integron_file && integron_file.size() > 0) ?
    "--integron_file ${integron_file}" : ""

  """
  set -euo pipefail

  # Get absolute paths for files
  # amr_path=\$(readlink -f "${amr_tsv}")
  # contigs_path=\$(readlink -f "${contigs_report}")

  python3 ${projectDir}/bin/compare_query.py \\
    --db_path ${db_in} \\
    --fasta_name ${meta.id} \\
    ${organism_value ? "--organism \"${organism_value}\"" : ""} \\
    --mash_dist_output "${mash_dist_output}" \\
    --amrfinder_output "${amr_tsv}" \\
    --contigs_report_path "${contigs_report}" \\
    ${ice_hits && ice_hits.size() > 0 ? "--filtered_hits_report_path ${ice_hits}" : ""} \\
    ${phage_coords && phage_coords.size() > 0 ? "--phage_report_path ${phage_coords}" : ""} \\
    ${compArg} \\
    ${tn3Arg} \\
    ${integronArg} \\
    --number ${params.number ?: 5} \\
    --output_format ${params.output_format ?: 'json'} \\
    --max_distance ${params.max_distance}
          
    2>&1 | tee register.log

    """
    stub:
    """
    touch ${meta.id}_resistome_differences.json
    touch analysis.log
    """
  }
