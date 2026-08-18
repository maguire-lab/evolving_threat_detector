process FILTER_DIAMOND_HITS {

    input:
    tuple val(meta), path(file)
    val tool_name
    val max_evalue
    val min_coverage

    output:
    tuple val(meta), path ("${meta.id}_${tool_name}_filtered_hits.tsv"), emit: filtered_iceberg_hits

    script:
    """
    awk -v ev="${max_evalue}" -v cov="${min_coverage}" 'BEGIN { OFS="\\t" }
        \$13 <= ev && \$4 > 0 && \$5 > 0 &&
	(\$6 / \$5) >= cov && (\$6 / \$4) >= cov { print }' "${file}" > "${meta.id}_${tool_name}_filtered_hits.tsv"
    """

   stub:
   """
   touch ${meta.id}_${tool_name}_filtered_hits.tsv
   """
}
