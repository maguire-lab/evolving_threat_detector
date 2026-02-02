process PARSE_GBK {
    
    tag "$meta.id"
    label 'process_low'
    
    container "quay.io/biocontainers/biopython:1.84"    

    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("*_contigs.fasta"), path("*_proteins.fasta"), path("*.gff3")

    script:
    """
    python ${projectDir}/bin/parse_gbk.py ${gbk} ${meta.id}
    """

    stub:
    """
    touch ${meta.id}_contigs.fasta
    touch ${meta.id}_proteins.fasta
    touch ${meta.id}.gff3
    """
}
