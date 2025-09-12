process INSERT_DB {
  tag "${meta.id}"
  label 'process_low'

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
  // build an optional GBK argument only if we actually have files
  def gbkArg = (gbk_files && gbk_files.size() > 0) ? 
    "--comp_gbk_files ${gbk_files.collect{ it.toString() }.join(' ')}" : ""
  """
  set -euo pipefail

  # Get absolute paths for files
  amr_path=\$(readlink -f "${amr_tsv}")
  contigs_path=\$(readlink -f "${contigs_report}")
  sketch_path=\$(readlink -f "${sketch_msh}")

  python3 ${projectDir}/bin/aggregate_output.py \\
    --db_path etd.db \\
    --fasta_name ${meta.id} \\
    --organism "${meta.organism}" \\
    --sketch_path "${sketch_msh}" \\
    --amrfinder_output "${amr_tsv}" \\
    --contigs_report_path "${contigs_report}" \\
    ${ice_hits && ice_hits.size() > 0 ? "--filtered_hits_report_path ${ice_hits}" : ""} \\
    ${phage_coords && phage_coords.size() > 0 ? "--phage_report_path ${phage_coords}" : ""} \\
    ${gbkArg} \\
    2>&1 | tee register.log

  # Capture any debug files created
  ls *_debug.txt 2>/dev/null || true
  """
}
