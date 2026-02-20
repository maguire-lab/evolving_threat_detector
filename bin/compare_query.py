#!/usr/bin/env python3

import sqlite3
import json
from pathlib import Path
import polars as pl
import argparse
from Bio import SeqIO
import sys
from aggregate_output import parse_plasmid_annotation, parse_ice_annotation, parse_phage_annotation, parse_composite_transposon_annotation,parse_tn3_transposon_annotation, parse_integron_annotation, is_proximal

def retrieve_closest_relatives(mash_dist_output, number=5):
    try:
        closest_relatives = []
        with open(mash_dist_output) as f:
            for line in f:
                fields = line.strip()
                if not line:
                    continue
                fields = line.split('\t')
                genome_name = fields[0].replace('_contigs.fasta','')
                mash_distance = float(fields[2])
                p_value = float(fields[3])
                numerator, denominator = map(int, fields[4].split('/'))
                similarity = numerator / denominator
                closest_relatives.append((genome_name, mash_distance, p_value, similarity))

        closest_relatives.sort(key=lambda x: x[1])
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
                start = int(parts[2])
                stop = int(parts[3])
                gene_name = parts[5]
                annotation = parts[6]
                query_amr_annotations.append({
                    "contig_id": contig_id,
                    "protein_id": protein_id,
                    "start": start,
                    "stop": stop,
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
                "contig_id": annotation["contig_id"],
                "start": annotation["start"],
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
                "contig_id": annotation["contig_id"],
                "start": annotation["start"],
                "ice_annotation": info["ice"],
                "ice_id": info["protein_id"]
                })

    return amr_ice_annotations

def merge_amr_phage_annotation(query_amr_annotations, phage_annotations, max_distance=5000):
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
        regions = phage_annotations.get(annotation["contig_id"])
        if regions:
            proximal_regions = []
            for r in regions:
                if is_proximal(annotation["start"], annotation["stop"],
                               r["start"], r["end"], max_distance):
                    proximal_regions.append(
                        f"{r['prophage_id']} ({r['start']}-{r['end']})"
                    )
            if proximal_regions:
                amr_phage_annotations.append({
                    "amr_gene": annotation["amr_gene"],
                    "contig_id": annotation["contig_id"],
                    "start": annotation["start"],
                    "prophage_annotation": ", ".join(proximal_regions)
                })

    return amr_phage_annotations

def merge_amr_comp_transposon(query_amr_annotations, composite_transposon_annotation, max_distance=5000):
    amr_comp_transposon_annotations = []
    for annotation in query_amr_annotations:
        info = composite_transposon_annotation.get(annotation["contig_id"])
        if info and info.get("start") is not None and info.get("end") is not None:
            if is_proximal(annotation["start"], annotation["stop"],
                           info["start"], info["end"], max_distance):
                amr_comp_transposon_annotations.append({
                    "contig_id": annotation["contig_id"],
                    "start": annotation["start"],
                    "amr_gene": annotation["amr_gene"],
                    "composite_transposon_annotation": info.get("composite_transposon_annotation")
                })
    return amr_comp_transposon_annotations

def merge_amr_tn3_transposon(query_amr_annotations, tn3_annotations):
    """Join AMR with Tn3 by contig_id"""
    amr_tn3_annotations = []
    for annotation in query_amr_annotations:
        info = tn3_annotations.get(annotation["contig_id"])
        if info:
            amr_tn3_annotations.append({
                "contig_id": annotation["contig_id"],
                "start": annotation["start"],
                "amr_gene": annotation["amr_gene"],
                "tn3_annotation": info.get("tn3_annotation")
            })
    return amr_tn3_annotations

def merge_amr_integron_annotation(query_amr_annotations, integron_annotations, max_distance=5000):
    amr_integron_annotations = []
    for annotation in query_amr_annotations:
        regions = integron_annotations.get(annotation["contig_id"])
        if regions:
            proximal_regions = []
            for r in regions:
                if is_proximal(annotation["start"], annotation["stop"],
                               r["start"], r["end"], max_distance):
                    proximal_regions.append(
                        f"{r['integron_id']} ({r['start']}-{r['end']})"
                    )
            if proximal_regions:
                amr_integron_annotations.append({
                    "amr_gene": annotation["amr_gene"],
                    "contig_id": annotation["contig_id"],
                    "start": annotation["start"],
                    "integron_annotation": ", ".join(proximal_regions)
                })
    return amr_integron_annotations

def merge_all_query_annotations(query_amr_annotations, plasmid_amr_annotations=None, phage_amr_annotations=None, ice_amr_annotations=None, comp_transposon_amr_annotations=None, tn3_amr_annotations=None, integron_amr_annotations=None):
      # Create the base DataFrame
      query_amr_annotations_df = pl.DataFrame(
          query_amr_annotations,
          schema=['contig_id', 'protein_id', 'start', 'stop', 'amr_gene', 'annotation'],
          orient="row"
      )

      # Start with the base DataFrame
      merged_df = query_amr_annotations_df

      # Join with plasmid annotations if provided
      if plasmid_amr_annotations:
          plasmid_df = pl.DataFrame(
              plasmid_amr_annotations,
              schema=['amr_gene', 'contig_id', 'start', 'plasmid_annotation', 'plasmid_id'],
              orient="row"
          )
          merged_df = merged_df.join(
                  plasmid_df.select(['amr_gene', 'contig_id', 'start', 'plasmid_annotation']), 
                  on=['amr_gene', 'contig_id', 'start'],
                  how='left'
                  )

      # Join with phage annotations if provided
      if phage_amr_annotations:
          phage_df = pl.DataFrame(
              phage_amr_annotations,
              schema=['amr_gene', 'contig_id', 'start', 'prophage_annotation'],
              orient="row"
          )
          merged_df = merged_df.join(
                  phage_df.select, 
                  on=['amr_gene', 'contig_id', 'start'], 
                  how='left'
                  )

      # Join with ICE annotations if provided
      if ice_amr_annotations:
          ice_df = pl.DataFrame(
              ice_amr_annotations,
              schema=['amr_gene', 'contig_id', 'start', 'ice_annotation', 'ice_id'],
              orient="row"
          )
          merged_df = merged_df.join(
                  ice_df.select(['amr_gene', 'contig_id', 'start', 'ice_annotation']),
                  on=['amr_gene', 'contig_id', 'start'],
                  how='left'
                  )

      # Join with composite transposon annotations if provided
      if comp_transposon_amr_annotations:
          comp_transposon_df = pl.DataFrame(
              comp_transposon_amr_annotations,
              schema=['contig_id', 'start', 'amr_gene', 'composite_transposon_annotation'],
              orient="row"
          )
          # Rename to avoid duplicate column names
          merged_df = merged_df.join(
              comp_transposon_df,
              on=['amr_gene', 'contig_id', 'start'],
              how='left',
              suffix='_comp'
          )

      # Join with Tn3 transposon annotations if provided
      if tn3_amr_annotations:
          tn3_df = pl.DataFrame(
                 tn3_amr_annotations,
                 schema=['contig_id','start', 'amr_gene', 'tn3_annotation'],
                 orient="row"
                 )
         # Rename to avoid duplicate column names
          merged_df = merged_df.join(
                 tn3_df,
                 on=['amr_gene', 'contig_id', 'start'],
                 how='left',
                 suffix='_tn3'
         )

    # Join with integron annotations if provided
      if integron_amr_annotations:
          integron_df = pl.DataFrame(
            integron_amr_annotations,
            schema=['amr_gene', 'contig_id', 'start', 'integron_annotation'],
            orient="row"
        )
          merged_df = merged_df.join(
            integron_df,
            on=['amr_gene', 'contig_id', 'start'],
            how='left'
        )
      else:
          merged_df = merged_df.with_columns(pl.lit(None).alias('integron_annotation'))


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
          (row['amr_gene'], row['contig_id'], row['start']): row for row in merged_df.to_dicts()
      }

      for genome_name, mash_distance, p_value, similarity in closest_relatives:
          cursor.execute(
              """
              SELECT g.genome_name, a.gene_name, a.contig_id, a.start, a.stop, a.amr_annotation, 
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
              schema=['genome_name', 'gene_name', 'contig_id', 'start', 'stop',
                     'amr_annotation', 'plasmid_annotation','integron_annotation', 
                     'prophage_annotation','composite_transposon_annotation',
                     'tn3_transposon_annotation', 'ice_annotation'],
              orient="row"
          )

          # Convert to dict for easier lookup
          closest_annotations_dict = {
              (row['gene_name'], row['contig_id'], row['start']): row for row in closest_annotations_df.to_dicts()
          }

          # Build gene-name sets for presence/absence comparison
          query_gene_names = set(key[0] for key in query_annotations_dict)
          relative_gene_names = set(key[0] for key in closest_annotations_dict)

          # Check for gained genes (genes in query but not in relative)
          gained_genes = query_gene_names - relative_gene_names
          for key, query_row in query_annotations_dict.items():
              gene = key[0]
              if gene in gained_genes:
                  differences.append({
                      'relative_genome': genome_name,
                      'gene_name': gene,
                      'mash_distance': mash_distance,
                      'p_value': p_value,
                      'similarity': similarity,

                      # Query annotations
                      'query_amr_annotation': query_row.get('annotation') or 'None',
                      'query_plasmid_annotation': query_row.get('plasmid_annotation') or 'None',
                      'query_prophage_annotation': query_row.get('prophage_annotation') or 'None',
                      'query_ice_annotation': query_row.get('ice_annotation') or 'None',
                      'query_composite_transposon_annotation':query_row.get('composite_transposon_annotation') or 'None',
                      'query_tn3_transposon_annotation': query_row.get('tn3_annotation') or 'None',

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
          lost_genes = relative_gene_names - query_gene_names
          for key, closest_row in closest_annotations_dict.items():
              gene = key[0]
              if gene in lost_genes:
                  differences.append({
                      'relative_genome': genome_name,
                      'gene_name': gene,
                      'mash_distance': mash_distance,
                      'p_value': p_value,
                      'similarity': similarity,

                      # Query annotations (all "Gene not present" since gene is absent)
                      'query_amr_annotation': 'Gene not present',
                      'query_plasmid_annotation': 'Gene not present',
                      'query_prophage_annotation': 'Gene not present',
                      'query_ice_annotation': 'Gene not present',
                      'query_composite_transposon_annotation': 'Gene not present',
                      'query_tn3_transposon_annotation': 'Gene not present',

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
    parser.add_argument('--comp_txt_files', type=Path, nargs='*', default=None, help='Paths to composite transposon TXT files')
    parser.add_argument('--tn3_txt_files', type=Path, nargs='*', default=None, help='Paths to Tn3 transposon TXT files')
    parser.add_argument('--integron_file', type=Path, default=None, help='Path to IntegronFinder .integrons file')
    parser.add_argument('--number', type=int, default=5, help='Number of closest genomes to consider. Default is 5.')
    parser.add_argument('--output_format', choices=['json', 'dataframe'], default='json', help='Output format for the resistome differences. Default is json.')
    parser.add_argument('--max_distance', type=int, default=5000, help='Maximum distance (bp) between AMR gene and MGE element to consider them co-located (default is 5000)')


    args = parser.parse_args()

    # Initialize variables
    plasmid_amr_annotations = None
    ice_amr_annotations = None
    phage_amr_annotations = None
    comp_transposon_amr_annotations = None
    tn3_amr_annotations = None
    integron_amr_annotations = None

    # Normalize comp_gbk_files to a list of strings (existing only)
    comp_txt_list = []
    if args.comp_txt_files:
        for p in args.comp_txt_files:
            p = Path(p)
            if p.exists():
                comp_txt_list.append(str(p))

    # Normalize tn3_gbk_files to a list of strings (existing only)
    tn3_txt_list = []
    if args.tn3_txt_files:
        for p in args.tn3_txt_files:
            p = Path(p)
            if p.exists():
                tn3_txt_list.append(str(p))


    db_path = Path(args.db_path)
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    closest_relatives =  retrieve_closest_relatives(args.mash_dist_output, args.number)

     # parse query amr_annotations
    query_amr_annotations = []
    if args.amrfinder_output and Path(args.amrfinder_output).exists():
        query_amr_annotations = parse_query_amr_annotation(args.amrfinder_output)
    else:
        print("Note: No AMR file provided or not found; skipping AMR insert.")

    # parse and merge plasmid annotations
    if query_amr_annotations and args.contigs_report_path and Path(args.contigs_report_path).exists():
        plasmid_annotations = parse_plasmid_annotation(args.contigs_report_path)
        plasmid_amr_annotations = merge_amr_plasmid_annotation(query_amr_annotations, plasmid_annotations)

    # parse and merge ICE annotations
    if query_amr_annotations and args.filtered_hits_report_path and Path(args.filtered_hits_report_path).exists():
        ice_annotations = parse_ice_annotation(args.filtered_hits_report_path)
        ice_amr_annotations = merge_amr_ice_annotation(query_amr_annotations, ice_annotations)

    # parse and merge prophage annotations
    if query_amr_annotations and args.phage_report_path and Path(args.phage_report_path).exists():
        phage_annotations = parse_phage_annotation(args.phage_report_path)
        phage_amr_annotations = merge_amr_phage_annotation(query_amr_annotations, phage_annotations, args.max_distance)

    # parse and merge composite transposon annotations
    if query_amr_annotations and comp_txt_list:
        comp_transposon_annotations = parse_composite_transposon_annotation(comp_txt_list)
        comp_transposon_amr_annotations = merge_amr_comp_transposon(query_amr_annotations, comp_transposon_annotations, args.max_distance)

    # parse and merge Tn3 transposon annotations
    if query_amr_annotations and tn3_txt_list:
        tn3_annotations = parse_tn3_transposon_annotation(tn3_txt_list)
        tn3_amr_annotations = merge_amr_tn3_transposon(query_amr_annotations, tn3_annotations)

    # parse and merge integron annotations
    if query_amr_annotations and args.integron_file and Path(args.integron_file).exists():
        integron_annotations = parse_integron_annotation(args.integron_file)
        if integron_annotations:
            integron_amr_annotations = merge_amr_integron_annotation(query_amr_annotations, integron_annotations, args.max_distance)
    
    merged_df =  merge_all_query_annotations(query_amr_annotations, plasmid_amr_annotations, phage_amr_annotations, ice_amr_annotations, comp_transposon_amr_annotations, tn3_amr_annotations, integron_amr_annotations)

    differences = compare_amr_annotations(cursor, merged_df, closest_relatives)

    prepare_output(differences, args.fasta_name, args.output_format)

    conn.close()


if __name__ == "__main__":
    main()
