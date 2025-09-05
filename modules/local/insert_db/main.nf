process INSERT_DB {
  tag "${meta.id}"
  label 'process_low'

  input:
  // per-genome bundle prepared by joins/mapping
  tuple val(meta), path(amr_tsv), path(contigs_report),
        path(ice_hits)      optional: true,
        path(phage_coords)  optional: true,
        path(gbk_files)     optional: true, collect: true
  
  // global (single) files:
  path sketch_msh
  path db_in

  output:
  path "etd.db"       , emit: db
  path "register.log" , emit: log

  script:
  // build an optional GBK argument only if we actually have files
  def gbkArg = (gbk_files && gbk_files.size() > 0) ? 
    "--comp-gbk-files ${gbk_files.collect{ it.toString() }.join(' ')}" : ""
  """
  set -euo pipefail

  # Seed or create DB
  if [ -s "${db_in}" ]; then cp "${db_in}" etd.db; else : > etd.db; fi

  # Collect GBK files if they exist
  # gbk_files=\$(ls *.gbk 2>/dev/null | tr '\\n' ' ' || echo "")


  bin/aggregate_output.py \\
    --db-path etd.db \\
    --fasta-name ${meta.id} \\
    --organism "${meta.organism}" \\
    --sketch-path "${sketch_msh}" \\
    --amrfinder-output "${amr_tsv}" \\
    --contigs-report-path "${contigs_report}" \\
    ${ ice_hits     ? "--filtered-hits-report-path ${ice_hits}" : ""} \\
    ${ phage_coords ? "--phage-report-path ${phage_coords}"     : ""} \\
    ${ (gbkArg } \\
    2>&1 | tee register.log
  """
}
