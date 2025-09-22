process TNCOMP_FINDER {
    tag "$meta.id"
    label 'medium'

    input:
    tuple val(meta), path(fasta)
    
    output:
    tuple val(meta), path("${fasta.baseName}_composite.txt"), optional: true, emit: report
    tuple val(meta), path("${fasta.baseName}_composite_*.gbk"), optional: true, emit: gbk
    path "blastn", emit: blastn
    path "info.txt", emit: info

    script:
    """
    python3 ${projectDir}/bin/tncomp/TnComp_finder.py \\
        -f $fasta \\
        -p ${task.cpus} \\
        -g -e 500 -k
    """

   stub:
   """
   touch ${fasta.baseName}_composite.txt
   touch ${fasta.baseName}_composite_candidate1.gbk
   touch ${fasta.baseName}_composite_candidate2.gbk
   touch blastn
   touch info.txt
   """
}
