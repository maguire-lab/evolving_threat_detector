import sqlite3
import json
from pathlib import Path
import polars as pl

def retrieve_closest_relatives(query_fasta, mash_dist_output, cursor, number=5):
    try:
        closest_relatives = []
        with open(mash_dist_output) as f:
            for line in f:
                fields = line.strip().split('\t')
                genome_name = Path(fields[0]).name
                p_value = float(fields[2])
                numerator, denominator = map(int, fields[4].split('/'))
                dist = numerator / denominator
                closest_genomes.append((genome_name, p_value, dist))

        closest_relatives.sort(key=lambda x: x[2])  # Sort by dist
        closest_relatives = closest_relatives[:number]
    except subprocess.CalledProcessError as e:
        print(f"Error running Mash dist: {e}")
        print(f"Command stderr: {e.stderr.decode()}")
        raise
    return closest_relatives

def compare_amr_annotations(cursor, query_annotations, closest_genomes):
    """
    Compare AMR annotations of the query genome against the closest genomes in the database.

    Parameters:
    cursor (sqlite3.Cursor): Database cursor object.
    query_annotations (list): Annotations of the query genome.
    closest_genomes (list): List of tuples (genome_name, p_value, dist) of closest genome names, distance metric, and their p-values.

    Returns:
    pl.DataFrame: A Polars DataFrame containing the resistome differences.
    """
    differences = []
    query_annotations_dict = {annotation['gene']: annotation['annotation'] for annotation in query_annotations if annotation['gene'] != 'Gene symbol'}

    for genome_name, p_value, dist in closest_genomes:
        cursor.execute(
            """
            SELECT g.genome_name, a.gene_name, a.amr_annotation, a.plasmid_annotation
            FROM annotations a
            JOIN genomes g ON a.genome_id = g.id
            WHERE g.genome_name = ? AND a.gene_name IS NOT NULL
            """,
            (genome_name,)
            )
        closest_annotations = cursor.fetchall()

        # Create a Polars DataFrame for the closest annotations
        closest_annotations_df = pl.DataFrame(
                closest_annotations,
                schema=['genome_name', 'gene_name', 'amr_annotation', 'plasmid_annotation'],
                orient="row"
                )

        # Compare each gene in the query genome
        for gene, query_annot in query_annotations_dict.items():
            # Filter closest annotations for the current gene
            gene_row = closest_annotations_df.filter(pl.col('gene_name') == gene)
            if gene_row.is_empty():
                closest_amr_annot = "Not present in relative"
                plasmid_annot = None

                # if query_annot != closest_annot:
                differences.append({
                    'fasta_name': genome_name,
                    'gene_name': gene,
                    'query_annotation': query_annot,
                    'closest_relatives_annotation': closest_amr_annot,
                    'difference_type': 'gained_gene',
                    'plasmid': plasmid_annot
                    })

        # Compare genes in closest genome to detect lost genes
        for row in closest_annotations_df.iter_rows(named=True):
            gene = row['gene_name']
            if gene not in query_annotations_dict:
                differences.append({
                    'fasta_name': genome_name,
                    'gene_name': gene,
                    'query_annotation': "Not present in query",
                    'closest_relatives_annotation': row['amr_annotation'],
                    'difference_type': 'lost_gene',
                    'plasmid': row['plasmid_annotation']
                    })

    return pl.DataFrame(differences)

