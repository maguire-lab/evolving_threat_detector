process GET_ICEBERG {
    label 'process_single'
    label 'error_retry_delay'

    output:
    path "ICE_aa_experimental_reformatted.fas", emit: iceberg

    script:
    """
    curl https://bioinfo-mml.sjtu.edu.cn/ICEberg2/download/ICE_aa_experimental.fas --output ICE_aa_experimental.fas

    sed -E 's/>(ICEberg\\|[0-9]+)\\s+(gi.*\\|)\\s+(.*)\\s+\\[(.*)\\]/>\\1_\\3_\\2_[\\4]/' <ICE_aa_experimental.fas | \
    tr ' ' '_' > ICE_aa_experimental_reformatted.fas
    """
    stub:
    """
    touch ICE_aa_experimental_reformatted.fas
    """
}
