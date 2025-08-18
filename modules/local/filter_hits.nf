process FILTER_DIAMOND_HITS {

    input:
    path blast_results
    val tool_name
    val min_pident
    val min_length

    output:
    path "${tool_name}_filtered_hits.tsv", emit: filtered_iceberg_hits

    script:
    """
    awk -v pid="${min_pident}" -v len="${min_length}" 'BEGIN { OFS="\\t" }
        \$3 >= pid && \$4 >= len { print }' "${blast_results}" > "${tool_name}_filtered_hits.tsv"
    """
}
