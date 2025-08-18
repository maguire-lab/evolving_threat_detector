process TN3_FINDER {
    tag "$meta.id"
    label 'medium'

    conda:
        "${projectDir}/modules/local/tn3finder/environment.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}.txt"), emit: report
    tuple val(meta), path("${fasta.baseName}_*.gbk"), optional: true, emit: gbk
    path "tblastn", emit: tblastn
    path "info.txt", emit: info

    script:
    """
    Tn3+TA_finder.py \\
        -f $fasta \\
        -g \\
        -t \${task.cpus}
    """
}
