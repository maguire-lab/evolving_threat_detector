process AMRFINDERPLUS_RUN {
    tag "$genomeID"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.12.8--h283d18e_0':
        'biocontainers/ncbi-amrfinderplus:3.12.8--h283d18e_0' }"

    input:
    tuple val(genomeID), path(genome)
    path db

    output:
    tuple val(genomeID), path("${prefix}.tsv")          , emit: report
    tuple val(genomeID), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    path "versions.yml"                             , emit: versions
    env VER                                         , emit: tool_version
    env DBVER                                       , emit: db_version

    when:
    task.ext.when == null || task.ext.when

    publishDir "${params.outdir}/amrfinderplus", mode: params.publish_dir_mode

    script:
    def args = task.ext.args ?: ''
    def is_compressed_fasta = genome.getName().endsWith(".gz") ? true : false
    def is_compressed_db = db.getName().endsWith(".gz") ? true : false
    prefix = task.ext.prefix ?: "${genomeID}"
    organism_param = params.organism ? "--organism ${params.organism} --mutation_all ${prefix}-mutations.tsv" : ""
    fasta_name = genome.getName().replace(".gz", "")
    fasta_param = "-n"
   // if (meta.containsKey("is_proteins")) {
       // if (meta.is_proteins) {
           // fasta_param = "-p"
       // }
   // }
    """
    if [ "$is_compressed_fasta" == "true" ]; then
        gzip -c -d $genome > $fasta_name
    fi

    if [ "$is_compressed_db" == "true" ]; then
        mkdir amrfinderdb
        tar xzvf $db -C amrfinderdb
    else
        mv $db amrfinderdb
    fi

    amrfinder \\
        $fasta_param $fasta_name \\
        $organism_param \\
        $args \\
        --database amrfinderdb \\
        --threads $task.cpus > ${prefix}.tsv

    VER=\$(amrfinder --version)
    DBVER=\$(echo \$(amrfinder --database amrfinderdb --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(echo \$(amrfinder --database amrfinderdb --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv

    VER=\$(amrfinder --version)
    DBVER=stub_version

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: stub_version
    END_VERSIONS
    """
}
