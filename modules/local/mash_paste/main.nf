process MASH_PASTE {
  label 'process_low'
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
      'biocontainers/mash:2.3--he348c14_1' }"

  // We take a *list* of files via .collect()
  input:
  path(msh_files) // this will be a list of *.msh from collect()

  output:
  path("all_genomes.msh"), emit: reference
  path("versions.yml")    , emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  # Write the file list to avoid cmdline length issues
  printf "%s\\n" ${msh_files} | sort -u > msh.list

  # mash paste <outprefix> <sketch1.msh> <sketch2.msh> ...
  # Use xargs to expand safely even for long lists
  xargs -a msh.list mash paste all_genomes

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      mash: \$(mash --version 2>&1)
  END_VERSIONS
  """

  stub:
  """
  touch all_genomes.msh
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      mash: stub
  END_VERSIONS
  """
}
