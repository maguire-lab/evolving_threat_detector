process PROKKA {
    tag "${genomeID}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3a/3af46b047c8fe84112adeaecf300878217c629b97f111f923ecf327656ddd141/data' :
        'community.wave.seqera.io/library/prokka_openjdk:10546cadeef11472' }"

    input:
    tuple val(genomeID), path(genome)
    //path proteins
    //path prodigal_tf

    output:
    tuple val(genomeID), path("${prefix}/*.gff"), emit: gff
    tuple val(genomeID), path("${prefix}/*.gbk"), emit: gbk
    tuple val(genomeID), path("${prefix}/*.fna"), emit: fna
    tuple val(genomeID), path("${prefix}/*.faa"), emit: faa
    tuple val(genomeID), path("${prefix}/*.ffn"), emit: ffn
    tuple val(genomeID), path("${prefix}/*.sqn"), emit: sqn
    tuple val(genomeID), path("${prefix}/*.fsa"), emit: fsa
    tuple val(genomeID), path("${prefix}/*.tbl"), emit: tbl
    tuple val(genomeID), path("${prefix}/*.err"), emit: err
    tuple val(genomeID), path("${prefix}/*.log"), emit: log
    tuple val(genomeID), path("${prefix}/*.txt"), emit: txt
    tuple val(genomeID), path("${prefix}/*.tsv"), emit: tsv
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args             = task.ext.args   ?: ''
    prefix               = task.ext.prefix ?: "${genomeID}"
    def input            = genome.toString() - ~/\.gz$/
    def decompress       = genome.getExtension() == "gz" ? "gunzip -c ${fasta} > ${input}" : ""
    def cleanup          = genome.getExtension() == "gz" ? "rm ${input}" : ""
    """
    ${decompress}

    prokka \\
        ${args} \\
        --cpus ${task.cpus} \\
        --prefix ${prefix} \\
        ${input}

    ${cleanup}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${genomeID}"
    """
    mkdir ${prefix}
    touch ${prefix}/${prefix}.gff
    touch ${prefix}/${prefix}.gbk
    touch ${prefix}/${prefix}.fna
    touch ${prefix}/${prefix}.faa
    touch ${prefix}/${prefix}.ffn
    touch ${prefix}/${prefix}.sqn
    touch ${prefix}/${prefix}.fsa
    touch ${prefix}/${prefix}.tbl
    touch ${prefix}/${prefix}.err
    touch ${prefix}/${prefix}.log
    touch ${prefix}/${prefix}.txt
    touch ${prefix}/${prefix}.tsv
    touch ${prefix}/${prefix}.gff

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        prokka: \$(echo \$(prokka --version 2>&1) | sed 's/^.*prokka //')
    END_VERSIONS
    """
}
