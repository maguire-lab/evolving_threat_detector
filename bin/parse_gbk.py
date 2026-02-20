#!/usr/bin/env python3
"""
Convert GenBank (.gbk) files to GFF3, contig nucleotide FASTA, and protein FASTA files.

This script uses BioPython to parse GenBank files and extract:
1. GFF3 annotations
2. Contig/chromosome nucleotide sequences
3. Protein sequences from CDS features

Requirements:
    pip install biopython

Usage:
    python parse_gbk.py input.gbk
"""

import sys
import os
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from datetime import datetime


def write_gff3_header(gff_file, source="GenBank"):
    """Write GFF3 format header with metadata."""
    gff_file.write("##gff-version 3\n")
    gff_file.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
    gff_file.write(f"##source-version {source}\n")


def format_gff3_attributes(feature):
    """
    Format GFF3 attributes from BioPython feature qualifiers.
    
    Chain of thought:
    1. Extract relevant qualifiers from the feature
    2. Format them according to GFF3 standards
    3. Handle special cases like gene names, products, etc.
    """
    attributes = []
    
    # Get feature ID - use locus_tag, gene, or generate one
    feature_id = None
    if 'locus_tag' in feature.qualifiers:
        feature_id = feature.qualifiers['locus_tag'][0]
    elif 'gene' in feature.qualifiers:
        feature_id = feature.qualifiers['gene'][0]
    else:
        # Generate ID based on location
        feature_id = f"{feature.type}_{feature.location.start}_{feature.location.end}"
    
    attributes.append(f"ID={feature_id}")
    
    # Add common attributes
    if 'gene' in feature.qualifiers:
        attributes.append(f"Name={feature.qualifiers['gene'][0]}")
    
    if 'product' in feature.qualifiers:
        # Escape special characters in product names
        product = feature.qualifiers['product'][0].replace(',', '%2C').replace(';', '%3B').replace('=', '%3D')
        attributes.append(f"product={product}")
    
    if 'note' in feature.qualifiers:
        note = ';'.join(feature.qualifiers['note']).replace(',', '%2C').replace(';', '%3B').replace('=', '%3D')
        attributes.append(f"Note={note}")
    
    return ';'.join(attributes)


def convert_strand(bio_strand):
    """Convert BioPython strand notation to GFF3 format."""
    if bio_strand == 1:
        return '+'
    elif bio_strand == -1:
        return '-'
    else:
        return '.'


def gbk_to_files(gbk_file, output_prefix=None):
    """
    Convert GenBank file to GFF3, nucleotide FASTA, and protein FASTA files.
    
    Chain of thought approach:
    1. Parse the GenBank file using BioPython
    2. For each record (contig/chromosome):
       a. Extract nucleotide sequence for FASTA
       b. Process features for GFF3 and protein extraction
    3. Write output files in appropriate formats
    
    Args:
        gbk_file (str): Path to input GenBank file
        output_prefix (str): Prefix for output files (optional)
    """
    
    # Determine output prefix
    if output_prefix is None:
        output_prefix = os.path.splitext(os.path.basename(gbk_file))[0]
    
    # Output file names
    gff_file = f"{output_prefix}.gff3"
    nucleotide_fasta = f"{output_prefix}_contigs.fasta"
    protein_fasta = f"{output_prefix}_proteins.fasta"
    
    print(f"Converting {gbk_file}...")
    print(f"Output files will be:")
    print(f"  - GFF3: {gff_file}")
    print(f"  - Nucleotide FASTA: {nucleotide_fasta}")
    print(f"  - Protein FASTA: {protein_fasta}")
    
    # Parse GenBank file and process each record
    with open(gff_file, 'w') as gff_out, \
         open(nucleotide_fasta, 'w') as nucl_out, \
         open(protein_fasta, 'w') as prot_out:
        
        # Write GFF3 header
        write_gff3_header(gff_out)
        
        record_count = 0
        feature_count = 0
        protein_count = 0
        
        for record in SeqIO.parse(gbk_file, "genbank"):
            record_count += 1
            
            # Write contig sequence directive to GFF3
            gff_out.write(f"##sequence-region {record.id} 1 {len(record.seq)}\n")
            
            # Write nucleotide sequence to FASTA
            nucl_out.write(f">{record.id}")
            if record.description:
                nucl_out.write(f" {record.description}")
            nucl_out.write(f"\n")
            
            # Write sequence in 80-character lines
            seq_str = str(record.seq)
            for i in range(0, len(seq_str), 80):
                nucl_out.write(seq_str[i:i+80] + "\n")
            
            # Process features for GFF3 and protein extraction
            for feature in record.features:
                if feature.type == "source":
                    continue

                feature_count += 1

                start = int(feature.location.start) + 1
                end = int(feature.location.end)
                strand = convert_strand(feature.location.strand)

                # Generate feature ID
                if feature.type == "CDS":
                    protein_count += 1
                    feature_id = f"{output_prefix}_{protein_count:05d}"
                elif 'locus_tag' in feature.qualifiers:
                    feature_id = feature.qualifiers['locus_tag'][0]
                elif 'gene' in feature.qualifiers:
                    feature_id = feature.qualifiers['gene'][0]
                else:
                    feature_id = f"{feature.type}_{start}_{end}"

                # Build GFF3 attributes
                attrs = [f"ID={feature_id}"]
                if 'gene' in feature.qualifiers:
                    attrs.append(f"Name={feature.qualifiers['gene'][0]}")
                if 'product' in feature.qualifiers:
                    product = feature.qualifiers['product'][0].replace(',', '%2C').replace(';', '%3B').replace('=', '%3D')
                    attrs.append(f"product={product}")

                gff_line = [
                    record.id, "GenBank", feature.type,
                    str(start), str(end), ".", strand, ".",
                    ";".join(attrs)
                ]
                gff_out.write("\t".join(gff_line) + "\n")

                # Extract protein sequence from CDS features
                if feature.type == "CDS" and 'translation' in feature.qualifiers:
                    prot_out.write(f">{feature_id}")
                    if 'product' in feature.qualifiers:
                        prot_out.write(f" {feature.qualifiers['product'][0]}")
                    prot_out.write(f" [location={record.id}:{start}..{end}({strand})]")
                    prot_out.write(f"\n")

                    prot_seq = feature.qualifiers['translation'][0]
                    for i in range(0, len(prot_seq), 80):
                        prot_out.write(prot_seq[i:i+80] + "\n")
    
    print(f"\nConversion completed successfully!")
    print(f"Processed {record_count} record(s)")
    print(f"Extracted {feature_count} feature(s)")
    print(f"Extracted {protein_count} protein sequence(s)")


def main():
    """Main function to handle command line arguments and file conversion."""
    if len(sys.argv) < 2:
        print("Usage: python gbk_converter.py <input.gbk> [output_prefix]")
        print("\nExample:")
        print("  python gbk_converter.py genome.gbk")
        print("  python gbk_converter.py genome.gbk my_genome")
        sys.exit(1)
    
    gbk_file = sys.argv[1]
    
    # Check if input file exists
    if not os.path.exists(gbk_file):
        print(f"Error: Input file '{gbk_file}' not found.")
        sys.exit(1)
    
    # Get output prefix if provided
    output_prefix = sys.argv[2] if len(sys.argv) > 2 else None
    
    try:
        gbk_to_files(gbk_file, output_prefix)
    except Exception as e:
        print(f"Error during conversion: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
