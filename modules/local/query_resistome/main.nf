process QUERY_RESISTOME {
  tag "${meta.id}"
  label 'process_medium'

  publishDir "${params.outdir}/resistome_analysis", mode: 'copy'

  input:
  // per-genome bundle prepared by joins/mapping
  tuple val(meta), path(amr_tsv), path(contigs_report),
        path(ice_hits),
        path(phage_coords),
        path(gbk_files),
        path(mash_dist_output)
  
  // global (single) files:
  path db_in

  output:
  path "${meta.id}_resistome_differences.*"       , emit: differences
  path "*.log"                , emit: logs, optional: true
  
  when:
  task.ext.when == null || task.ext.when


  script:
  // build an optional GBK argument only if we actually have files
  def gbkArg = (gbk_files && gbk_files.size() > 0) ? 
    "--comp_gbk_files ${gbk_files.collect{ it.toString() }.join(' ')}" : ""
  """
  set -euo pipefail

  # Get absolute paths for files
  # amr_path=\$(readlink -f "${amr_tsv}")
  # contigs_path=\$(readlink -f "${contigs_report}")

  python3 ${projectDir}/bin/compare_query.py \\
    --db_path ${db_in} \\
    --fasta_name ${meta.id} \\
    --organism "${meta.organism}" \\
    --mash_dist_output "${mash_dist_output}" \\
    --amrfinder_output "${amr_tsv}" \\
    --contigs_report_path "${contigs_report}" \\
    ${ice_hits && ice_hits.size() > 0 ? "--filtered_hits_report_path ${ice_hits}" : ""} \\
    ${phage_coords && phage_coords.size() > 0 ? "--phage_report_path ${phage_coords}" : ""} \\
    ${gbkArg} \\
    --number ${params.number ?: 5} \\
    --output_format ${params.output_format ?: 'json'} \\
          
    2>&1 | tee register.log

    """
    stub:
    """
    touch ${meta.id}_resistome_differences.json
    touch analysis.log
    """
  }
