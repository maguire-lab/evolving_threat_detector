process MOBSUITE_RECON {
    tag "$genomeID"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mob_suite:3.1.9--pyhdfd78af_0':
        'biocontainers/mob_suite:3.1.9--pyhdfd78af_0' }"

    input:
    tuple val(genomeID), path(genome)

    output:
    tuple val(genomeID), path("results/$genomeID/chromosome.fasta")    , emit: chromosome
    tuple val(genomeID), path("results/$genomeID/contig_report.txt")   , emit: contig_report
    tuple val(genomeID), path("results/$genomeID/plasmid_*.fasta")     , emit: plasmids        , optional: true
    tuple val(genomeID), path("results/$genomeID/mobtyper_results.txt"), emit: mobtyper_results, optional: true
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when


    publishDir "${params.outdir}/mob_output", mode: params.publish_dir_mode

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${genomeID}"
    def is_compressed = genome.getName().endsWith(".gz") ? true : false
    def fasta_name = genome.getName().replace(".gz", "")
    """
    # rerunning mob recon
    if [ "$is_compressed" == "true" ]; then
        gzip -c -d $genome > $fasta_name
    fi

    mob_recon \\
        --infile $fasta_name \\
        $args \\
        --num_threads $task.cpus \\
        --outdir "results/$prefix" \\
        --sample_id $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p results

    touch results/chromosome.fasta
    touch results/contig_report.txt
    touch results/plasmids

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mobsuite: \$(echo \$(mob_recon --version 2>&1) | sed 's/^.*mob_recon //; s/ .*\$//')
    END_VERSIONS
    """
}
