process PARSE_GBK {
    
    tag "$meta.id"
    label 'process_low'
    conda = './environment.yml'
    
    input:
    tuple val(meta), path(gbk)

    output:
    tuple val(meta), path("*_contigs.fasta"), path("*_proteins.fasta"), path("*.gff3")

    script:
    """
    python3 ${projectDir}/bin/parse_gbk.py ${gbk} ${meta.id}
    """
}
