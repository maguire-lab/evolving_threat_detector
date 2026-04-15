process TN3_FINDER {
    tag "$meta.id"
    label 'process_high'

    container 'quay.io/biocontainers/tn3_ta_finder:1.0.1--hdfd78af_0'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${fasta.baseName}.txt"), optional: true, emit: report
    tuple val(meta), path("${fasta.baseName}_*.gbk"), optional: true, emit: gbk
    path "tblastn", emit: tblastn
    path "info.txt", emit: info

    script:
    """
    Tn3+TA_finder.py \\
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
