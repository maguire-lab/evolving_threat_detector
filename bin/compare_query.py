#!/usr/bin/env python3

import sqlite3
import json
from pathlib import Path
import polars as pl
import argparse
from Bio import SeqIO
import sys
from aggregate_output import parse_plasmid_annotation, parse_ice_annotation, parse_phage_annotation, parse_composite_transposon_annotation

def retrieve_closest_relatives(mash_dist_output, number=5):
    try:
        closest_relatives = []
        with open(mash_dist_output) as f:
            for line in f:
                fields = line.strip().split('\t')
                #genome_name = Path(fields[0]).name
                genome_name = fields[0].split('_')[0]
                p_value = float(fields[2])
                numerator, denominator = map(int, fields[4].split('/'))
                dist = 1 - (numerator / denominator)
                closest_relatives.append((genome_name, p_value, dist))

        closest_relatives.sort(key=lambda x: x[2])
        closest_relatives = closest_relatives[:number]
    except Exception as e:
        print(f"Error reading Mash file: {e}")
        raise
    return closest_relatives

def parse_query_amr_annotation(amrfinder_output):
    """
    Read the AMRFinder output (TSV), pick columns by index, and return a list of dicts.

    Parameters:
    amrfinder_output (str): Path to the output file from upstream amrfinder module.
    genome_id (int): The ID of the genome entry.

    Returns:
    amr_annotations (list): A list of annotations where each annotation is a dictionary with keys 'genome_id', 'gene' and 'annotation'.
    """
    try:
        query_amr_annotations = []
        with open(amrfinder_output, "r") as f:
            # skip header
            next(f, None)
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 7:
                    # skip malformed/short lines
                    continue
                contig_id = parts[1].split()[0]
                protein_id = parts[0]
                gene_name = parts[5]
                annotation = parts[6]
                query_amr_annotations.append({
                    "contig_id": contig_id,
                    "protein_id": protein_id,
                    "amr_gene": gene_name,
                    "annotation": annotation
                })
        return query_amr_annotations
    except Exception as e:
        print(f"An unexpected error occurred while parsing {amrfinder_output}: {e}")
        raise

def merge_amr_plasmid_annotation(query_amr_annotations, plasmid_annotations):
    """
    Join AMR rows with plasmid annotations by contig_id.

    Returns:
        list of dicts: {
            "genome_id": ...,
            "amr_gene": ...,
            "plasmid_annotation": ...,
        }
    """
    amr_plasmid_annotations = []
    for annotation in query_amr_annotations:
        info = plasmid_annotations.get(annotation["contig_id"])
        if info:
            amr_plasmid_annotations.append({
                "amr_gene": annotation["amr_gene"],
                "plasmid_annotation": info["plasmid"],
                "plasmid_id": info["contig_id"]
                })

    return amr_plasmid_annotations

def merge_amr_ice_annotation(query_amr_annotations, ice_annotations):
    """
    Join AMR rows with ice annotations by contig_id.

    Returns:
        list of dicts: {
            "genome_id": ...,
            "amr_gene": ...,
            "ice_annotation": ...,
        }
    """
    amr_ice_annotations = []
    for annotation in query_amr_annotations:
        info = ice_annotations.get(annotation["protein_id"])
        if info:
            amr_ice_annotations.append({
                "amr_gene": annotation["amr_gene"],
                "ice_annotation": info["ice"],
                "ice_id": info["protein_id"]
                })

    print(f"Created {len(amr_ice_annotations)} AMR-ICE annotations")

    return amr_ice_annotations

def merge_amr_phage_annotation(query_amr_annotations, phage_annotations):
    """
    Join AMR rows with prophage annotations by contig_id.

    Returns:
        list of dicts: {
            "genome_id": ...,
            "amr_gene": ...,
            "prophage_annotation": ...,
        }
    """
    amr_phage_annotations = []
    for annotation in query_amr_annotations:
        info = phage_annotations.get(annotation["contig_id"])
        if info:
            amr_phage_annotations.append({
                "amr_gene": annotation["amr_gene"],
                "prophage_annotation": info["prophage"],
                "phage_id": info["contig_id"]
                })

    return amr_phage_annotations

def merge_amr_comp_transposon(query_amr_annotations, composite_transposon_annotation):
    amr_comp_transposon_annotations = []
    for annotation in query_amr_annotations:
        info = composite_transposon_annotation.get(annotation["contig_id"])
        if info:
            amr_comp_transposon_annotations.append({
                "contig_id": annotation["contig_id"],
                "amr_gene":  annotation["amr_gene"],
                "composite_transposon_annotation": info.get("composite_transposon_annotation")
            })
    return amr_comp_transposon_annotations

def merge_all_query_annotations(query_amr_annotations, plasmid_amr_annotations=None, phage_amr_annotations=None, ice_amr_annotations=None, comp_transposon_amr_annotations=None):
      # Create the base DataFrame
      query_amr_annotations_df = pl.DataFrame(
          query_amr_annotations,
          schema=['contig_id', 'protein_id', 'amr_gene', 'annotation']
      )

      # Start with the base DataFrame
      merged_df = query_amr_annotations_df

      # Join with plasmid annotations if provided
      if plasmid_amr_annotations:
          plasmid_df = pl.DataFrame(
              plasmid_amr_annotations,
              schema=['plasmid_id', 'amr_gene', 'plasmid_annotation']
          )
          merged_df = merged_df.join(plasmid_df, on='amr_gene', how='left')

      # Join with phage annotations if provided
      if phage_amr_annotations:
          phage_df = pl.DataFrame(
              phage_amr_annotations,
              schema=['phage_id', 'amr_gene', 'prophage_annotation']
          )
          merged_df = merged_df.join(phage_df, on='amr_gene', how='left')

      # Join with ICE annotations if provided
      if ice_amr_annotations:
          ice_df = pl.DataFrame(
              ice_amr_annotations,
              schema=['ice_id', 'amr_gene', 'ice_annotation']
          )
          merged_df = merged_df.join(ice_df, on='amr_gene', how='left')

      # Join with composite transposon annotations if provided
      if comp_transposon_amr_annotations:
          comp_transposon_df = pl.DataFrame(
              comp_transposon_amr_annotations,
              schema=['contig_id', 'amr_gene', 'composite_transposon_annotation']
          )
          # Rename to avoid duplicate column names
          merged_df = merged_df.join(
              comp_transposon_df.rename({'contig_id': 'transposon_contig_id'}),
              on='amr_gene',
              how='left'
          )

      return merged_df

def compare_amr_annotations(cursor, merged_df, closest_relatives):
      """
      Compare AMR gene presence/absence between query genome and closest relatives.

      Parameters:
      cursor (sqlite3.Cursor): Database cursor object.
      merged_df (pl.DataFrame): Merged annotations of the query genome.
      closest_relatives (list): List of tuples (genome_name, p_value, dist) of closest genome 
  names, distance metric, and their p-values.

      Returns:
      pl.DataFrame: A Polars DataFrame containing the resistome differences.
      """
      differences = []

      # Convert merged_df to dictionary for easier lookup
      query_annotations_dict = {
          row['amr_gene']: row for row in merged_df.to_dicts()
      }

      for genome_name, p_value, dist in closest_relatives:
          cursor.execute(
              """
              SELECT g.genome_name, a.gene_name, a.amr_annotation, 
                     a.plasmid_annotation, a.integron_annotation, a.prophage_annotation, 
                     a.composite_transposon_annotation, a.tn3_transposon_annotation, 
  a.ice_annotation
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
              schema=['genome_name', 'gene_name', 'amr_annotation', 'plasmid_annotation',
                     'integron_annotation', 'prophage_annotation',
  'composite_transposon_annotation',
                     'tn3_transposon_annotation', 'ice_annotation']
          )

          # Convert to dict for easier lookup
          closest_annotations_dict = {
              row['gene_name']: row for row in closest_annotations_df.to_dicts()
          }

          # Check for gained genes (genes in query but not in relative)
          for gene, query_row in query_annotations_dict.items():
              if gene not in closest_annotations_dict:
                  differences.append({
                      'relative_genome': genome_name,
                      'gene_name': gene,
                      'p_value': p_value,
                      'distance': dist,

                      # Query annotations
                      'query_amr_annotation': query_row.get('annotation') or 'None',
                      'query_plasmid_annotation': query_row.get('plasmid_annotation') or 'None',
                      'query_prophage_annotation': query_row.get('prophage_annotation') or 'None',
                      'query_ice_annotation': query_row.get('ice_annotation') or 'None',
                      'query_composite_transposon_annotation':
  query_row.get('composite_transposon_annotation') or 'None',

                      # Relative annotations (all "Gene not present" since gene is absent)
                      'relative_amr_annotation': 'Gene not present',
                      'relative_plasmid_annotation': 'Gene not present',
                      'relative_integron_annotation': 'Gene not present',
                      'relative_prophage_annotation': 'Gene not present',
                      'relative_composite_transposon_annotation': 'Gene not present',
                      'relative_tn3_transposon_annotation': 'Gene not present',
                      'relative_ice_annotation': 'Gene not present',

                      'difference_type': 'gained_gene'
                  })

          # Check for lost genes (genes in relative but not in query)
          for gene, closest_row in closest_annotations_dict.items():
              if gene not in query_annotations_dict:
                  differences.append({
                      'relative_genome': genome_name,
                      'gene_name': gene,
                      'p_value': p_value,
                      'distance': dist,

                      # Query annotations (all "Gene not present" since gene is absent)
                      'query_amr_annotation': 'Gene not present',
                      'query_plasmid_annotation': 'Gene not present',
                      'query_prophage_annotation': 'Gene not present',
                      'query_ice_annotation': 'Gene not present',
                      'query_composite_transposon_annotation': 'Gene not present',

                      # Relative annotations
                      'relative_amr_annotation': closest_row.get('amr_annotation') or 'None',
                      'relative_plasmid_annotation': closest_row.get('plasmid_annotation') or 'None',
                      'relative_integron_annotation': closest_row.get('integron_annotation') or
  'None',
                      'relative_prophage_annotation': closest_row.get('prophage_annotation') or
  'None',
                      'relative_composite_transposon_annotation':
  closest_row.get('composite_transposon_annotation') or 'None',
                      'relative_tn3_transposon_annotation':
  closest_row.get('tn3_transposon_annotation') or 'None',
                      'relative_ice_annotation': closest_row.get('ice_annotation') or 'None',

                      'difference_type': 'lost_gene'
                  })

      return pl.DataFrame(differences)

def prepare_output(differences, fasta_name, output_format='json'):
    """
    Function to report resistome differences between a query genome and its closest reference sequences.

    Parameters:
    differences: A Polars DataFrame containing the resistome differences.
    output_format (str): Output format for the result ('json' or 'dataframe'). Default is 'json'.

    Returns:
    str or pl.DataFrame: A JSON string or a Polars DataFrame containing the resistome differences.
    """
    if output_format.lower() == 'json':
        filename = f"{fasta_name}_resistome_differences.json"
        differences_dict = differences.to_dicts()
        with open(filename, 'w') as f:
            json.dump(differences_dict, f, indent=2)
        return filename
    elif output_format.lower() == 'dataframe':
        filename = f"{fasta_name}_resistome_differences.csv"
        differences.write_csv(filename)
        return filename
    else:
        raise ValueError(f"Unsupported output format: {output_format}. Use 'json' or 'dataframe'.")


def main():
    parser = argparse.ArgumentParser(description='Queries etd db for closest relatives and parses query genome annotations.')
    parser.add_argument('--db_path', type=Path, required=True, help='Path to the existing SQLite database.')
    parser.add_argument('--fasta_name', type=str, required=True, help='Genome fasta id')
    parser.add_argument('--organism', type=str, default=None, help='Organism name')
    parser.add_argument('--mash_dist_output', type=Path, help='Path to the mash dist output file')
    parser.add_argument('--amrfinder_output', type=Path, help='Path to amrfinderplus TSV file')
    parser.add_argument('--contigs_report_path', type=Path, help='Path to mobsuite contigs_report.txt')
    parser.add_argument('--filtered_hits_report_path', type=Path, default=None, help='Path to to ICE filtered_hits TSV')
    parser.add_argument('--phage_report_path', type=Path, default=None, help='Path to the prophage report TSV')
    parser.add_argument('--comp_gbk_files', type=Path, nargs='*', default=None, help='Paths to composite transposon GBK files (one per candidate)')
    parser.add_argument('--number', type=int, default=5, help='Number of closest genomes to consider. Default is 5.')
    parser.add_argument('--output_format', choices=['json', 'dataframe'], default='json', help='Output format for the resistome differences. Default is json.')


    args = parser.parse_args()

    # Initialize variables
    plasmid_amr_annotations = None
    ice_amr_annotations = None
    phage_amr_annotations = None
    comp_transposon_amr_annotations = None

    # Normalize comp_gbk_files to a list of strings (existing only)
    gbk_list = []
    if args.comp_gbk_files:
        for p in args.comp_gbk_files:
            p = Path(p)
            if p.exists():
                gbk_list.append(str(p))


    db_path = Path(args.db_path)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    closest_relatives =  retrieve_closest_relatives(args.mash_dist_output, args.number)

     # parse and store amr_annotations
    query_amr_annotations = []
    if args.amrfinder_output and Path(args.amrfinder_output).exists():
        query_amr_annotations = parse_query_amr_annotation(args.amrfinder_output)
    else:
        print("Note: No AMR file provided or not found; skipping AMR insert.")

    # parse and store plasmid annotations
    if query_amr_annotations and args.contigs_report_path and Path(args.contigs_report_path).exists():
        plasmid_annotations = parse_plasmid_annotation(args.contigs_report_path)
        plasmid_amr_annotations = merge_amr_plasmid_annotation(query_amr_annotations, plasmid_annotations)

    # parse and store ICE annotations
    if query_amr_annotations and args.filtered_hits_report_path and Path(args.filtered_hits_report_path).exists():
        ice_annotations = parse_ice_annotation(args.filtered_hits_report_path)
        ice_amr_annotations = merge_amr_ice_annotation(query_amr_annotations, ice_annotations)

    # parse and store prophage annotations
    if query_amr_annotations and args.phage_report_path and Path(args.phage_report_path).exists():
        phage_annotations = parse_phage_annotation(args.phage_report_path)
        phage_amr_annotations = merge_amr_phage_annotation(query_amr_annotations, phage_annotations)

    # parse and store composite transposon annotations
    if query_amr_annotations and gbk_list:
        comp_transposon_annotations = parse_composite_transposon_annotation(gbk_list)
        comp_transposon_amr_annotations = merge_amr_comp_transposon(query_amr_annotations, comp_transposon_annotations)
    
    merged_df =  merge_all_query_annotations(query_amr_annotations, plasmid_amr_annotations, phage_amr_annotations, ice_amr_annotations, comp_transposon_amr_annotations)

    differences = compare_amr_annotations(cursor, merged_df, closest_relatives)

    prepare_output(differences, args.fasta_name, args.output_format)

    conn.close()


if __name__ == "__main__":
    main()
