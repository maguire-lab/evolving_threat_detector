process FILTER_DIAMOND_HITS {

    input:
    tuple val(meta), path(file)
    val tool_name
    val min_pident
    val min_length

    output:
    tuple val(meta), path ("${tool_name}_filtered_hits.tsv"), emit: filtered_iceberg_hits

    script:
    """
    awk -v pid="${min_pident}" -v len="${min_length}" 'BEGIN { OFS="\\t" }
        \$3 >= pid && \$6 >= len { print }' "${file}" > "${tool_name}_filtered_hits.tsv"
    """
}
