process MASH_PASTE_BATCH {
  label 'process_low'
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
      'biocontainers/mash:2.3--he348c14_1' }"

  input:
  path(msh_files)

  output:
  path("batch_*.msh"), emit: batch_msh
  path("versions.yml"), emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  # Write all sketch file paths to a list (one per line, deduplicated)
  ls *.msh | sort -u > msh.list

  # Paste this chunk into a single intermediate sketch file
  # task.index gives each parallel batch a unique name
  # The -l flag reads input paths from the list file in one invocation,
  # avoiding xargs splitting and the "file already exists" error
  mash paste batch_\${RANDOM} -l msh.list

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      mash: \$(mash --version 2>&1)
  END_VERSIONS
  """

  stub:
  """
  touch batch_0.msh
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
      mash: stub
  END_VERSIONS
  """
}

process MASH_PASTE_FINAL {
  label 'process_low'
  conda "${moduleDir}/environment.yml"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
      'biocontainers/mash:2.3--he348c14_1' }"

  input:
  path(batch_files)

  output:
  path("all_genomes.msh"), emit: reference
  path("versions.yml"),    emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  """
  # Write all intermediate batch sketch paths to a list
  ls batch_*.msh | sort -u > msh.list

  # Paste all batch intermediates into the final combined sketch
  mash paste all_genomes -l msh.list

  # Clean up intermediate batch files
  rm batch_*.msh

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
