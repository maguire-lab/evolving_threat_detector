// NOTE: Custom version of MASH_SKETCH — includes publishDir directive and custome input and o//utput variable names

process MASH_SKETCH {
    tag "$genomeID"
    label 'process_medium'
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mash:2.3--he348c14_1' :
        'biocontainers/mash:2.3--he348c14_1' }"

    input:
    tuple val(genomeID), path(genome)

    output:
    tuple val(genomeID), path("*.msh")        , emit: mash
    tuple val(genomeID), path("*.mash_stats") , emit: stats
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    publishDir "${params.outdir}/mash", mode: params.publish_dir_mode

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${genomeID}"
    """
    mash \\
        sketch \\
        $args \\
        $genome \\
        -p $task.cpus \\
        -o ${prefix} \\
        2> >(tee ${prefix}.mash_stats >&2)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.msh
    touch ${prefix}.mash_stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mash: \$(mash --version 2>&1)
    END_VERSIONS
    """
}
