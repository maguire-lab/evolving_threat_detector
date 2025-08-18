process TNCOMP_FINDER {
    tag "$meta.id"
    label 'medium'

    conda:
        "${projectDir}/modules/local/tncomp_finder/environment.yml"

    input:
    tuple val(meta), path(fasta)
    
    output:
    tuple val(meta), path("${fasta.baseName}_composite.txt"), emit: report
    tuple val(meta), path("${fasta.baseName}_composite_*.gbk"), optional: true, emit: gbk
    path "blastn", emit: blastn
    path "info.txt", emit: info

    script:
    """
    python3 TnComp_finder.py \\
        -f $fasta \\
        -o . \\
        -p ${task.cpus} \\
        -g -e 500 -k
    """
}
