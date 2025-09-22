process TN3_FINDER {
    tag "$meta.id"
    label 'medium'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}.txt"), optional: true, emit: report
    tuple val(meta), path("${fasta.baseName}_*.gbk"), optional: true, emit: gbk
    path "tblastn", emit: tblastn
    path "info.txt", emit: info

    script:
    """
    python3 ${projectDir}/bin/tn3/Tn3+TA_finder.py \\
        -f $fasta \\
        -g \\
        -t ${task.cpus}
    """

    stub:
    """
    touch ${fasta.baseName}.txt
    touch ${fasta.baseName}_candidate1.gbk
    touch ${fasta.baseName}_candidate2.gbk
    touch tblastn
    touch info.txt
    """
}
