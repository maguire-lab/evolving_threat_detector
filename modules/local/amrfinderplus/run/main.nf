process AMRFINDERPLUS_RUN {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:4.2.7--hf69ffd2_0':
        'biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0' }"

    input:
    tuple val(meta), path(fasta), path(protein), path(gff)
    path db

    output:
    tuple val(meta), path("${prefix}.tsv")          , emit: report
    tuple val(meta), path("${prefix}-mutations.tsv"), emit: mutation_report, optional: true
    path "versions.yml"                             , emit: versions
    env VER                                         , emit: tool_version
    env DBVER                                       , emit: db_version

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def is_compressed_fasta   = fasta.getName().endsWith(".gz") ? true : false
    def is_compressed_protein = protein.getName().endsWith(".gz") ? true : false
    def is_compressed_gff     = gff.getName().endsWith(".gz") ? true : false  
    //def is_compressed_db      = db.getName().endsWith(".gz") ? true : false
    prefix = task.ext.prefix ?: "${meta.id}"
    organism_param = meta.organism ? "--organism ${meta.organism} --mutation_all ${prefix}-mutations.tsv" : ""
    fasta_name   = fasta.getName().replace(".gz", "")
    protein_name = protein.getName().replace(".gz", "")
    gff_name     = gff.getName().replace(".gz", "")
    annotation_format = gff_name.endsWith(".gff") ? "prokka" : "bakta"
    
   // fasta_param = "-n"
    //if (meta.containsKey("is_proteins")) {
        //if (meta.is_proteins) {
            //fasta_param = "-p"
        //}
    //}

    """
    if [ "$is_compressed_fasta" == "true" ]; then
        gzip -c -d $fasta > $fasta_name
    fi

    if [ "$is_compressed_protein" == "true" ]; then
        gzip -c -d $protein > $protein_name
    fi

    if [ "$is_compressed_gff" == "true" ]; then
        gzip -c -d $gff > $gff_name
    fi

    # use db already an extracted directory
    # ln -s ${db} amrfinderdb

    # combined amrfinderplus run

    amrfinder \\
        -n $fasta_name \\
        -p $protein_name \\
        --gff $gff_name \\
        --annotation_format $annotation_format \\
        $organism_param \\
        $args \\
        --database ${db}/latest \\
        --threads $task.cpus > ${prefix}.tsv

    VER=\$(amrfinder --version)
    DBVER=\$(echo \$(amrfinder --database ${db}/latest --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev)

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        amrfinderplus: \$(amrfinder --version)
        amrfinderplus-database: \$(echo \$(echo \$(amrfinder --database ${db}/latest --database_version 2> stdout) | rev | cut -f 1 -d ' ' | rev))
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
