#!/usr/bin/env python3

import sqlite3
from pathlib import Path
import argparse
from Bio import SeqIO
import sys
import json

DATABASE_PATH = 'etd.db'

def connect(db_path = DATABASE_PATH):
    """
    Open a SQLite connection with foreign keys enforced.
    Caller owns commit/close.
    """
    conn = sqlite3.connect(str(db_path))
    conn.execute("PRAGMA foreign_keys=ON")
    return conn

def init_db(db_path=DATABASE_PATH):
    """
    Initialize the database schema (create tables if they don't exist).
    
    Parameters:
    db_path (str): Path to the SQLite database file.
    """
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()

        # Create Genomes table
        cursor.execute('''CREATE TABLE IF NOT EXISTS genomes (
                          id INTEGER PRIMARY KEY,
                          genome_name TEXT NOT NULL UNIQUE,
                          organism TEXT)''')

        # Create Sketch table for storing the single sketch file
        cursor.execute('''CREATE TABLE IF NOT EXISTS sketch (
                          id INTEGER PRIMARY KEY,
                          sketch_path TEXT NOT NULL UNIQUE)''')

        # Create Annotations table
        cursor.execute('''CREATE TABLE IF NOT EXISTS annotations (
                          id INTEGER PRIMARY KEY,
                          genome_id INTEGER REFERENCES genomes(id),
                          gene_name TEXT,
                          contig_id TEXT,
                          start INTEGER,
                          stop INTEGER,
                          amr_annotation TEXT,
                          amr_output_path TEXT,
                          plasmid_annotation TEXT,
                          plasmid_output_path TEXT,
                          integron_annotation TEXT,
                          integron_output_path TEXT,
                          prophage_annotation TEXT,
                          prophage_output_path TEXT,
                          composite_transposon_annotation TEXT,
                          composite_transposon_output_path TEXT,
                          tn3_transposon_annotation TEXT,
                          tn3_transposon_output_path TEXT,
                          ice_annotation TEXT,
                          ice_output_path TEXT,
                          UNIQUE(genome_id, gene_name, contig_id, start))''')

        # Create MGE elements table — stores ALL detected MGEs per genome,
        # independent of whether they overlap with an AMR gene
        cursor.execute('''CREATE TABLE IF NOT EXISTS mge_elements (
                          id INTEGER PRIMARY KEY,
                          genome_id INTEGER REFERENCES genomes(id),
                          mge_type TEXT NOT NULL,
                          mge_name TEXT,
                          contig_id TEXT,
                          start_pos INTEGER,
                          end_pos INTEGER,
                          output_path TEXT,
                          UNIQUE(genome_id, mge_type, contig_id, start_pos, end_pos))''')

        # Create Contigs table — one row per contig per genome.
        # contig_length supports contig-edge (Indeterminate) detection;
        # is_circular supports wrap-around distance on closed replicons
        cursor.execute('''CREATE TABLE IF NOT EXISTS contigs (
                          id INTEGER PRIMARY KEY,
                          genome_id INTEGER REFERENCES genomes(id),
                          contig_id TEXT,
                          contig_length INTEGER,
                          is_circular INTEGER DEFAULT 0,
                          UNIQUE (genome_id, contig_id))''')


        # Index to speed up UPDATE/SELECT by (genome_id, gene_name)
        cursor.execute("""
        CREATE INDEX IF NOT EXISTS ix_annotations_gid_gene
        ON annotations(genome_id, gene_name, contig_id, start)
        """)

        # Index for speed up lookup by genome and MGE type
        cursor.execute("""
        CREATE INDEX IF NOT EXISTS ix_mge_elements_gid_type
        ON mge_elements(genome_id, mge_type)
        """)

        conn.commit()
        conn.close()
    except sqlite3.Error as e:
        print(f"Error initializing database: {e}")
        raise

def create_genome_entry(cursor, fasta_name, organism):
    """
    Create and store a Genome entry in the database if it does not exist, including the organism name.

    Parameters:
    cursor (sqlite3.Cursor): Database cursor object.
    fasta_name (str): The name of the FASTA file.
    organism (str): The name of the organism.

    Returns:
    int: The ID of the created genome entry.
    """
    try:
        # Try to find an existing row by (genome_name, organism)
        cursor.execute(
            "SELECT id FROM genomes WHERE genome_name = ? AND IFNULL(organism,'') = IFNULL(?, '')",
            (fasta_name, organism)
        )
        row = cursor.fetchone()
        if row:
            return int(row[0])
        # Else insert
        cursor.execute(
                "INSERT INTO genomes (genome_name, organism) VALUES (?, ?)", (fasta_name, organism))
        genome_id = cursor.lastrowid
        if genome_id is None:
            raise sqlite3.IntegrityError(f"Failed to insert genome '{fasta_name}'.")
        return genome_id
    except sqlite3.Error as e:
        print(f"Error inserting genome '{fasta_name}': {e}")
        raise

def store_sketch(sketch_path_published, cursor):
    """
    Store Mash sketch in the DB.

    Parameters:
    sketch_path (Path) : Path to the sketch file generated by upstream mash paste module
    cursor (sqlite3.Cursor): Database cursor object.
    """
    #sketch_path = Path(sketch_path_published).resolve()

    #if not sketch_path.exists():
       # raise FileNotFoundError(f"Sketch not found: {sketch_path}")

    #if sketch_path_published.suffix != ".msh":
       # print(f"Warning: expected .msh, got {sketch_path.suffix}: {sketch_path}")

    cursor.execute(
        "INSERT OR IGNORE INTO sketch (sketch_path) VALUES (?) ",
        (str(sketch_path_published),)
    )

def store_mge_elements(mge_records, cursor, output_path):
    """
    Insert all detected MGE elements into the mge_elements table.

    Parameters:
        mge_records (list): List of dicts, each with keys:
            genome_id, mge_type, mge_name, contig_id,
            start_pos, end_pos
        cursor: SQLite cursor
        output_path (str): Published path to the tool's output file
    """
    out_str = str(output_path) if output_path else None
    inserted = 0
    for rec in mge_records:
        try:
            cursor.execute(
                """
                INSERT OR IGNORE INTO mge_elements
                    (genome_id, mge_type, mge_name, contig_id,
                     start_pos, end_pos, output_path)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    rec["genome_id"],
                    rec["mge_type"],
                    rec.get("mge_name", None),
                    rec.get("contig_id", None),
                    rec.get("start_pos", None),
                    rec.get("end_pos", None),
                    out_str
                )
            )
            inserted += cursor.rowcount
        except sqlite3.IntegrityError:
            pass
    print(f"  MGE elements inserted: {inserted} ({rec.get('mge_type', 'unknown')})")
    return inserted


def parse_amr_annotation(amrfinder_output, genome_id):
    """
    Read the AMRFinder output (TSV), pick columns by index, and return a list of dicts.

    Parameters:
    amrfinder_output (str): Path to the output file from upstream amrfinder module.
    genome_id (int): The ID of the genome entry.

    Returns:
    amr_annotations (list): A list of annotations where each annotation is a dictionary with keys 'genome_id', 'gene' and 'annotation'.
    """
    try:
        amr_annotations = []
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
                amr_annotations.append({
                    "genome_id": genome_id,
                    "contig_id": contig_id,
                    "protein_id": protein_id,
                    "start": start,
                    "stop": stop,
                    "amr_gene": gene_name,
                    "annotation": annotation
                })
        return amr_annotations
    except Exception as e:
        print(f"An unexpected error occurred while parsing {amrfinder_output}: {e}")
        raise


def store_amr_annotation(amr_annotations, cursor, output_path):
    """
    Store AMR annotations in the database.

    Parameters:
    amr_annotations (list): list of amr_annotation dicts.
    output_path (str): Path to the output file from upstream amrfinder module.
    cursor (sqlite3.Cursor): Database cursor object..
    """
    
    out_str = str(output_path)
    try:
        for annotation in amr_annotations:
            cursor.execute(
                """
                INSERT INTO annotations 
                (genome_id, gene_name, contig_id, start, stop, amr_annotation, amr_output_path) 
                VALUES (?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(genome_id, gene_name, contig_id, start)
                DO UPDATE SET amr_annotation = excluded.amr_annotation,
                              amr_output_path = excluded.amr_output_path
                """,
                (
                    annotation['genome_id'], 
                    annotation['amr_gene'], 
                    annotation['contig_id'], 
                    annotation['start'], 
                    annotation['stop'], 
                    annotation['annotation'], 
                    out_str
                    )
                )
    except sqlite3.Error as e:
        print(f"Error storing AMR annotations: {e}")
        raise

def is_proximal(amr_start, amr_stop, element_start, element_end, max_distance, 
                contig_length=None, is_circular=False):
    """
    Test whether an AMR gene is within max_distance bp of a genomic element.

    Parameters:
        amr_start, amr_stop      (int): AMR gene coordinates, 1-based.
        element_start, element_end (int): element coordinates, 1-based.
        max_distance             (int): proximity window, bp.
        contig_length      (int, opt.): length of the shared contig. Required
                                        for wrap-around; ignored if absent.
        is_circular       (bool, opt.): True only for closed replicons. Must be
                                        False for draft contigs, where wrapping
                                        would be meaningless.

    Returns:
        bool: True if:
      - The AMR gene overlaps the element, OR
      - The gap between the nearest edges is <= max_distance

    All coordinates are 1-based genomic positions on the same contig.
    On a circular replicon the shorter of the two ways round is used.
    If contig_length is inconsistent with the coordinates, the wrap is
    discarded and only the linear gap is considered.
    """
    # Normalise so start <= stop (handles reverse strand)
    a_lo, a_hi = min(amr_start, amr_stop), max(amr_start, amr_stop)
    e_lo, e_hi = min(element_start, element_end), max(element_start, element_end)

    # If they overlap, distance is 0
    if a_lo <= e_hi and e_lo <= a_hi:
        return True

    # Otherwise, compute gap between nearest edges
    if a_lo > e_hi:                      # AMR gene lies to the right
        linear_gap = a_lo - e_hi
        wrap_gap = (contig_length - a_hi) + e_lo if contig_length else None
    else:                                # AMR gene lies to the left
        linear_gap = e_lo - a_hi
        wrap_gap = (contig_length - e_hi) + a_lo if contig_length else None

    if wrap_gap is not None and wrap_gap < 0:
        wrap_gap = None

    gap = (min(linear_gap, wrap_gap)
           if (is_circular and wrap_gap is not None) else linear_gap)
    return gap <= max_distance

def window_is_truncated(gene_start, gene_stop, max_distance,
                        contig_length, is_circular=False):
    """
    True when the proximity window extends past a contig end, so part of it
    lies in sequence the assembly does not contain. A circle has no ends.
    """
    if is_circular or not contig_length:
        return False
    lo, hi = min(gene_start, gene_stop), max(gene_start, gene_stop)
    return (lo - max_distance) < 1 or (hi + max_distance) > contig_length

def build_protein_contig_map(gbk_path, genome_id):
    """
    Parse a GenBank file to build a mapping from protein locus_tag
    to its contig location and coordinates.

    Parameters:
        gbk_path (str/Path): Path to the genome's GBK/GBFF file.
        genome_id (str): Genome identifier (e.g. "GCA_000290555.2"),
                         used as prefix for sequential protein IDs.

    Returns:
        dict: {protein_id: {
            "contig_id": str,  # the contig/record this CDS sits on
            "start": int,      # 1-based start position
            "end": int         # 1-based end position
        }}
    """
    protein_map = {}
    protein_count = 0
    for record in SeqIO.parse(str(gbk_path), 'genbank'):
        contig_id = record.id                   # e.g. "NZ_CP012345.1"
        for feat in record.features:
            if feat.type != "CDS":
                continue
            # Increment counter for every CDS (same order as parse_gbk.py)
            protein_count += 1
            # Build the same ID that parse_gbk.py writes to the protein FASTA
            prot_id = f"{genome_id}_{protein_count:05d}"
            protein_map[prot_id] = {
                "contig_id": contig_id,
                "start": int(feat.location.start) + 1,  # BioPython is 0-based
                "end": int(feat.location.end)
            }
    print(f'Built protein-contig map: {len(protein_map)} proteins from {gbk_path}')
    return protein_map

def build_contig_map(gbk_path):
    """
    Parse a GenBank file and return contig length and topology.

    Returns:
        dict: {contig_id: {"length": int, "circular": bool}}
    """
    contig_map = {}
    for record in SeqIO.parse(str(gbk_path), 'genbank'):
        topology = str(record.annotations.get('topology', 'linear')).lower()
        contig_map[record.id] = {
            "length": len(record.seq),
            "circular": topology == 'circular',
        }
    print(f'Built contig map: {len(contig_map)} contigs from {gbk_path}')
    return contig_map


def store_contigs(contig_map, genome_id, cursor):
    """Write contig lengths and topology for one genome."""
    rows = [(genome_id, cid, info["length"], 1 if info["circular"] else 0)
            for cid, info in contig_map.items()]
    cursor.executemany(
        "INSERT OR REPLACE INTO contigs "
        "(genome_id, contig_id, contig_length, is_circular) VALUES (?, ?, ?, ?)",
        rows)
    return len(rows)

def build_ice_element_metadata(iceberg_fasta_path):
    """
    Parse the raw ICEberg FASTA to extract element_id and functional
    category for each protein. The reformatted header format is:
    >ICEberg|{element_id}_{protein_name}_gi|{gi}|{db}|{accession}|_[{organism}]

    Parameters:
        iceberg_fasta_path (str/Path): Path to ICE_aa_experimental_reformatted.fas

    Returns:
        dict: {sseqid_string: {
            "element_id": str,     # e.g. "68"
            "protein_name": str,   # e.g. "int(Tn916)"
            "category": str        # one of: integrase, relaxase, t4cp, t4ss, cargo
        }}
    """
    # Keywords validated against ICEberg FASTA protein names
    INTEGRASE_KW = ['integrase', 'recombinase', 'xerc', 'xerd', 'excisionase']

    RELAXASE_KW  = ['relaxase', 'moba', 'mobb', 'mobc',
                    'mobilization_protein', 'mobilisation_protein', 'trai']

    T4CP_KW      = ['virb4', 'coupling', 'trae', 'vird4', 'trad']

    T4SS_KW      = ['virb', 'sex_pilus', 'pilus_assembly', 'mating_pair',
                    'conjugative_transfer', 'conjugal_transfer',
                    'type_iv_secret', 'type-iv_secret',
                    'type_iv_b_pilus', 'type_iv_pilus', 'type_4_pilus',
                    'traa', 'trab', 'traf', 'trah',
                    'trak', 'tral', 'trau', 'traw', 'traq',
                    'trbb', 'trbc', 'trbd', 'trbe', 'trbf', 'trbg',
                    'trbi', 'trbj', 'trbl',
                    'conjugation_signal_peptidase']

    metadata = {}
    for record in SeqIO.parse(str(iceberg_fasta_path), 'fasta'):
        header = record.description          # full FASTA header line
        sid = record.id                       # the part before first whitespace

        # Parse element_id from header: ICEberg|{element_id}_...
        element_id = 'unknown'
        protein_name = header
        try:
            # After reformatting: >ICEberg|68_int(Tn916)_gi|...
            after_pipe = header.split("|", 1)[1]   # "68_int(Tn916)_gi|..."
            element_id = after_pipe.split("_", 1)[0]  # "68"
            rest = after_pipe.split("_", 1)[1] if "_" in after_pipe else ""
            # protein_name is everything up to the next _gi|
            if "_gi|" in rest:
                protein_name = rest.split("_gi|")[0]
            else:
                protein_name = rest
        except (IndexError, ValueError):
            pass  # keep defaults

        # Classify into functional category based on protein name
        name_lower = protein_name.lower()
        category = "cargo"  # default: not a signature protein
        if any(kw in name_lower for kw in INTEGRASE_KW):
            category = "integrase"
        elif any(kw in name_lower for kw in RELAXASE_KW):
            category = "relaxase"
        elif any(kw in name_lower for kw in T4CP_KW):
            category = "t4cp"
        elif any(kw in name_lower for kw in T4SS_KW):
            category = "t4ss"
        # Special case: TraC — catch without matching tetracycline/extracellular/bacitracin
        elif 'trac' in name_lower and not any(x in name_lower for x in
                ['tetrac', 'extrac', 'intrac', 'bacitrac']):
            category = "t4ss"
        # Special case: TraN — only match at word boundaries to avoid transcription/transposase
        elif ('_tran_' in name_lower or name_lower.endswith('_tran') or
              name_lower.startswith('tran_') or name_lower == 'tran'):
            category = "t4ss"

        metadata[sid] = {
            "element_id": element_id,
            "protein_name": protein_name,
            "category": category
        }

    print(f'Built ICEberg metadata: {len(metadata)} proteins, '
          f'{len(set(m["element_id"] for m in metadata.values()))} unique elements')
    return metadata


def parse_plasmid_annotation(contigs_report_path):
    """
    Read in the mob-suite contigs_report.txt and return a dict of plasmid annotations

    Parameters:
    contigs_report_path (str): Path to the output file from upstream mob-suite module

    Returns:
    plasmid annotations (dict): dict of plamid annotations.
    """
    try:
        plasmid_annotations = {}
        with open(contigs_report_path) as f:
            next(f, None)
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue

                molecule_type_full = parts[1].strip()
                if not molecule_type_full.lower().startswith("plasmid"):
                # first_token = molecule_type_full.split(None, 1)[0].lower() if molecule_type_full else ""
                # if first_token != "plasmid":
                    continue

                plasmid = parts[2].strip()

                contig_id_norm = parts[4].split()[0]
                if not contig_id_norm:
                    continue

                plasmid_annotations[contig_id_norm] = {
                    "contig_id": contig_id_norm,
                    "plasmid": plasmid
                }

        if not plasmid_annotations:
            print(f"Warning: No plasmid contigs found in {contigs_report_path}")

        return plasmid_annotations

    except Exception as e:
        print(f"An unexpected error occurred while parsing {contigs_report_path}: {e}")
        raise

def merge_amr_plasmid_annotation(amr_annotations, plasmid_annotations):
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
    for annotation in amr_annotations:
        info = plasmid_annotations.get(annotation["contig_id"])
        if info:
            amr_plasmid_annotations.append({
                "genome_id": annotation["genome_id"],
                "contig_id": annotation["contig_id"],
                "start": annotation["start"],
                "amr_gene": annotation["amr_gene"],
                "plasmid_annotation": info["plasmid"],
                "plasmid_id": info["contig_id"]
                })

    return amr_plasmid_annotations

def store_amr_plasmid_annotations(amr_plasmid_annotations, cursor, output_path):
    """
    Update existing AMR rows with plasmid info.

    amr_plasmid_annotations: list of dicts like:
      {"genome_id": int, "contig_id": str, "amr_gene": str, "plasmid_annotation": str}

    Returns:
      int: number of table rows updated.
    """
    out_str = str(output_path)

    total_updated = 0
    for row in amr_plasmid_annotations:
        genome_id = row.get("genome_id")
        gene      = row.get("amr_gene")
        contig_id = row.get("contig_id")
        start     = row.get("start")
        plasmid   = row.get("plasmid_annotation")

        cursor.execute(
            """
            UPDATE annotations
                SET plasmid_annotation = ?,
                    plasmid_output_path = ?
            WHERE genome_id = ?
                AND gene_name = ?
                AND contig_id = ?
                AND start = ?
            """,
                (plasmid, out_str, genome_id, gene, contig_id, start)
            )
        total_updated += cursor.rowcount
    return total_updated

#def parse_ice_annotation(filtered_hits_report_path):
#    """
#    Read in the ICE filtered_hits.tsv report and return a dict of ICE annotations
#
#    Parameters:
#    filtered_hits_report_path (str): Path to the output file from upstream ICE module
#
#    Returns:
#    ICE annotations (dict): dict of ICE annotations.
#    """
#    try:
#        ice_annotations = {}
#        with open(filtered_hits_report_path) as f:
#            #next(f, None)
#            for line in f:
#                parts = line.strip().split("\t")
#                if len(parts) < 5:
#                    continue
#
#                protein_id = parts[0].strip()
#                if not protein_id:
#                    continue
#                annotation = parts[1].strip()
#
#                ice_annotations[protein_id] = {
#                    "protein_id": protein_id,
#                    "ice": annotation
#                }
#
#        if not ice_annotations:
#            print(f"Warning: No ice annotations found in {filtered_hits_report_path}")
#
#        return ice_annotations
#
#   except Exception as e:
#       print(f"An unexpected error occurred while parsing {filtered_hits_report_path}: {e}")
#        raise

def parse_ice_annotation(filtered_hits_report_path, protein_contig_map, ice_element_metadata):
    """
    Parse DIAMOND-vs-ICEberg filtered hits and produce element-level ICE
    annotations grouped by contig, suitable for proximity-based AMR matching.

    Parameters:
        filtered_hits_report_path: Path to DIAMOND filtered_hits.tsv
            Columns: qseqid  sseqid  pident  slen  qlen  length ...
        protein_contig_map: dict from build_protein_contig_map()
        ice_element_metadata: dict from build_ice_element_metadata()
    
    Returns:
        dict: {contig_id: [{
            "ice_label": str,     # e.g. "ICE_68 (integrase,relaxase,t4ss)"
            "element_id": str,    # e.g. "68"
            "start": int,         # leftmost coordinate of element hits on this contig
            "end": int,           # rightmost coordinate of element hits on this contig
            "categories": set     # e.g. {"integrase", "relaxase", "t4ss"}
        }, ...]}
    """
    try:
        # Parse hits and group by (contig_id, element_id)
        from collections import defaultdict
        groups = defaultdict(lambda: {'starts': [], 'ends': [], 'categories': set()})

        with open(filtered_hits_report_path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue

                qseqid = parts[0].strip()   # genome protein locus_tag
                sseqid = parts[1].strip()   # ICEberg protein ID

                # Look up which contig this query protein sits on
                prot_info = protein_contig_map.get(qseqid)
                if not prot_info:
                    continue  # protein not found in GBK

                # Look up the ICEberg element metadata for the subject
                ice_info = ice_element_metadata.get(sseqid)
                if not ice_info:
                    continue  # subject not in our metadata

                contig_id  = prot_info['contig_id']
                element_id = ice_info['element_id']
                category   = ice_info['category']
                key = (contig_id, element_id)
                groups[key]['starts'].append(prot_info['start'])
                groups[key]['ends'].append(prot_info['end'])
                groups[key]['categories'].add(category)

        # Filter groups — require >= 2 distinct SIGNATURE categories
        SIGNATURE_CATS = {'integrase', 'relaxase', 't4cp', 't4ss'}
        ice_annotations = defaultdict(list)

        for (contig_id, element_id), data in groups.items():
            sig_cats = data['categories'] & SIGNATURE_CATS
            if len(sig_cats) < 2:
                continue  # insufficient evidence for a real ICE element

            cat_str = ','.join(sorted(sig_cats))
            ice_label = f"ICE_{element_id} ({cat_str})"

            ice_annotations[contig_id].append({
                "ice_label": ice_label,
                "element_id": element_id,
                "start": min(data["starts"]),
                "end": max(data["ends"]),
                "categories": sig_cats
            })

        if not ice_annotations:
            print(f'Warning: No ICE elements passed the 2-category filter in {filtered_hits_report_path}')
        else:
            total = sum(len(v) for v in ice_annotations.values())
            print(f'Found {total} ICE elements across {len(ice_annotations)} contigs')

        return dict(ice_annotations)

    except Exception as e:
        print(f'An unexpected error occurred while parsing {filtered_hits_report_path}: {e}')
        raise



#def merge_amr_ice_annotation(amr_annotations, ice_annotations):
#    """
#    Join AMR rows with ice annotations by contig_id.

#    Returns:
#        list of dicts: {
#            "genome_id": ...,
#            "amr_gene": ...,
#            "ice_annotation": ...,
#        }
#    """
#    amr_ice_annotations = []
#    for annotation in amr_annotations:
#        info = ice_annotations.get(annotation["protein_id"])
#        if info:
#            amr_ice_annotations.append({
#                "genome_id": annotation["genome_id"],
#                "contig_id": annotation["contig_id"],
#                "start": annotation["start"],
#                "amr_gene": annotation["amr_gene"],
#                "ice_annotation": info["ice"],
#                "ice_id": info["protein_id"]
#                })
#
#    print(f"Created {len(amr_ice_annotations)} AMR-ICE annotations")
#
#    return amr_ice_annotations

def merge_amr_ice_annotation(amr_annotations, ice_annotations, max_distance=5000, contig_map=None):
    """
    Join AMR rows with ICE annotations by contig proximity.
    For each AMR gene, check all ICE element regions on the same contig.
    If any element boundary is within max_distance bp, annotate the AMR gene.

    Parameters:
        amr_annotations: list of dicts from parse_amr (with contig_id, start, stop)
        ice_annotations: dict from parse_ice_annotation (contig_id -> [regions])
        max_distance: int, bp threshold for proximity (default 5000)
        contig_map: dict from build_contig_map

    Returns:
        list of dicts: {genome_id, contig_id, start, amr_gene, ice_annotation}
    """
    amr_ice_annotations = []
    for annotation in amr_annotations:
        # Look up ICE elements on this AMR gene's contig
        regions = ice_annotations.get(annotation["contig_id"])
        if regions:
            cinfo = (contig_map or {}).get(annotation["contig_id"], {})
            proximal_elements = []
            for r in regions:
                # same is_proximal() function from above
                if is_proximal(annotation["start"], annotation["stop"],
                               r["start"], r["end"], max_distance,
                               cinfo.get("length"),             
                               cinfo.get("circular", False)):
                    proximal_elements.append(
                        f"{r['ice_label']} ({r['start']}-{r['end']})"
                    )
            if proximal_elements:
                amr_ice_annotations.append({
                    "genome_id": annotation["genome_id"],
                    "contig_id": annotation["contig_id"],
                    "start": annotation["start"],
                    "amr_gene": annotation["amr_gene"],
                    "ice_annotation": ", ".join(proximal_elements)
                })

    print(f'Created {len(amr_ice_annotations)} AMR-ICE annotations')
    return amr_ice_annotations


def store_amr_ice_annotations(amr_ice_annotations, cursor, output_path):
    """
    Update existing AMR rows with ice info.

    amr_ice_annotations: list of dicts like:
      {"genome_id": int, "contig_id": str, "amr_gene": str, "ice_annotation": str}

    Returns:
      int: number of table rows updated.
    """
    out_str = str(output_path)
    total_updated = 0
    for row in amr_ice_annotations:
        genome_id = row.get("genome_id")
        gene      = row.get("amr_gene")
        contig_id = row.get("contig_id")
        start     = row.get("start")
        ice   = row.get("ice_annotation")

        cursor.execute(
            """
            UPDATE annotations
                SET ice_annotation = ?,
                    ice_output_path = ?
            WHERE genome_id = ?
                AND gene_name = ?
                AND contig_id = ?
                AND start = ?
            """,
                (ice, str(out_str), genome_id, gene, contig_id, start)
            )
        total_updated += cursor.rowcount
    return total_updated

def parse_phage_annotation(phage_report_path):
    """
    Read in the PhiSpy prophage_coordinates.tsv and return a dict of prophage annotations.

    Each contig maps to a list of prophage regions with their boundaries.

    Parameters:
    phage_report_path (str): Path to the output file from upstream PhiSpy module

    Returns:
    dict: {contig_id: [{"prophage_id": str, "start": int, "end": int}, ...]}
    """
    try:
        phage_annotations = {}
        with open(phage_report_path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                parts = line.split("\t")
                if len(parts) < 4:
                    continue

                prophage_id = parts[0]
                contig_id = parts[1]
                if not contig_id:
                    continue
                phage_start = int(parts[2])
                phage_end = int(parts[3])

                if contig_id not in phage_annotations:
                    phage_annotations[contig_id] = []

                phage_annotations[contig_id].append({
                    "prophage_id": prophage_id,
                    "start": phage_start,
                    "end": phage_end
                })

        if not phage_annotations:
            print(f"Warning: No phage annotation found in {phage_report_path}")

        return phage_annotations

    except Exception as e:
        print(f"An unexpected error occurred while parsing {phage_report_path}: {e}")
        raise

def merge_amr_phage_annotation(amr_annotations, phage_annotations,
                               max_distance=5000, contig_map=None):
    """
    Join AMR rows with prophage annotations by contig_id.
    A contig may contain multiple prophage regions.

    Returns:
        list of dicts with genome_id, contig_id, start, amr_gene, prophage_annotation
    """
    amr_phage_annotations = []
    for annotation in amr_annotations:
        regions = phage_annotations.get(annotation["contig_id"])
        if regions:
            cinfo = (contig_map or {}).get(annotation["contig_id"], {})
            # Build annotation string from proximal prophage regions on this contig
            proximal_regions = []
            for r in regions:
                if is_proximal(annotation["start"], annotation["stop"],
                               r["start"], r["end"], max_distance,
                               cinfo.get("length"), 
                               cinfo.get("circular", False)):
                    proximal_regions.append(
                        f"{r['prophage_id']} ({r['start']}-{r['end']})"
                    )
            if proximal_regions:
                amr_phage_annotations.append({
                    "genome_id": annotation["genome_id"],
                    "contig_id": annotation["contig_id"],
                    "start": annotation["start"],
                    "amr_gene": annotation["amr_gene"],
                    "prophage_annotation": ", ".join(proximal_regions)
                })

    return amr_phage_annotations

def store_amr_phage_annotations(amr_phage_annotations, cursor, output_path):
    """
    Update existing AMR rows with prophage info.

    amr_phage_annotations: list of dicts like:
      {"genome_id": int, "contig_id": str, "amr_gene": str, "prophage_annotation": str}

    Returns:
      int: number of table rows updated.
    """
    out_str = str(output_path)
    total_updated = 0
    for row in amr_phage_annotations:
        genome_id = row.get("genome_id")
        gene      = row.get("amr_gene")
        contig_id = row.get("contig_id")
        start     = row.get("start")
        prophage   = row.get("prophage_annotation")

        cursor.execute(
            """
            UPDATE annotations
                SET prophage_annotation  = ?,
                    prophage_output_path = ?
            WHERE genome_id   = ?
                AND gene_name = ?
                AND contig_id = ?
                AND start     = ?
            """,
                (prophage, str(out_str), genome_id, gene, contig_id, start)
            )
        total_updated += cursor.rowcount
    return total_updated

def parse_composite_transposon_annotation(comp_txt_files):
    """
    Parse composite transposon TXT output files.

    Each file has:
      Line 1: QUERY: <contig_id> <organism> ...
      Then *Candidate N* blocks with a Feature table.
      Feature table columns: Feature, Position, Strand, Length(bp), Type

    Extract the IS element names (Type column) from rows where
    Feature contains 'transposon'.

    Returns:
        dict keyed by contig_id, each value is a dict with
        'contig_id' and 'composite_transposon_annotation' (comma-separated string).
    """
    composite_transposon_annotations = {}

    for txt_file in comp_txt_files:
        try:
            contig_id = None
            is_elements = []
            positions = []
            in_feature_table = False

            with open(txt_file) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        in_feature_table = False
                        continue

                    # Extract contig_id from QUERY line
                    if line.startswith("QUERY:"):
                        contig_id = line.split()[1]
                        continue

                    # Skip candidate headers, length/distance/order lines
                    if line.startswith("*Candidate") or line.startswith("Length:") or line.startswith("Distance:") or line.startswith("Order:"):
                        continue

                    # Detect the Feature table header
                    if line.startswith("Feature") and "Position" in line:
                        in_feature_table = True
                        continue

                    # Parse Feature table rows
                    if in_feature_table:
                        parts = line.split()
                        if len(parts) >= 5 and "transposon" in parts[1].lower():
                            # parts: [A/B, transposon, start..end, strand, length, type]
                            # But "A transposon" is split as ['A', 'transposon', ...]
                            is_type = parts[5] if len(parts) >= 6 else parts[4]
                            is_elements.append(is_type)

                            # Extract genomic positions
                            pos_str = parts[2]  # e.g. "14294..15082"
                            if ".." in pos_str:
                                pos_parts = pos_str.split("..")
                                positions.append(int(pos_parts[0]))
                                positions.append(int(pos_parts[1]))

                    # Stop parsing at sequence sections
                    if line.startswith(">") or line.startswith("Predicted ORFs"):
                        in_feature_table = False
                        break

            if contig_id and is_elements:
                if contig_id not in composite_transposon_annotations:
                    composite_transposon_annotations[contig_id] = {
                        "contig_id": contig_id,
                        "composite_transposon_annotation": [],
                        "start": None,
                        "end": None
                    }
                composite_transposon_annotations[contig_id]["composite_transposon_annotation"].extend(is_elements)
                # Update boundaries
                if positions:
                    existing_start = composite_transposon_annotations[contig_id]["start"]
                    existing_end = composite_transposon_annotations[contig_id]["end"]
                    new_start = min(positions)
                    new_end = max(positions)
                    composite_transposon_annotations[contig_id]["start"] = min(new_start, existing_start) if existing_start is not None else new_start
                    composite_transposon_annotations[contig_id]["end"] = max(new_end, existing_end) if existing_end is not None else new_end

        except Exception as e:
            print(f"Error parsing {txt_file}: {e}")

    # Deduplicate and convert to comma-separated strings
    for contig_id in composite_transposon_annotations:
        elements = composite_transposon_annotations[contig_id]["composite_transposon_annotation"]
        unique_elements = list(set(elements))
        composite_transposon_annotations[contig_id]["composite_transposon_annotation"] = ", ".join(unique_elements)

    return composite_transposon_annotations

def merge_amr_comp_transposon(amr_annotations, composite_transposon_annotation,
                              max_distance=5000, contig_map=None):
    amr_comp_transposon_annotations = []
    for annotation in amr_annotations:
        info = composite_transposon_annotation.get(annotation["contig_id"])
        if info and info.get("start") is not None and info.get("end") is not None:
            cinfo = (contig_map or {}).get(annotation["contig_id"], {})
            if is_proximal(annotation["start"], annotation["stop"],
                           info["start"], info["end"], max_distance,
                           cinfo.get("length"),
                           cinfo.get("circular", False)):
                amr_comp_transposon_annotations.append({
                    "genome_id": annotation["genome_id"],
                    "contig_id": annotation["contig_id"],
                    "start": annotation["start"],
                    "amr_gene": annotation["amr_gene"],
                    "composite_transposon_annotation": info.get("composite_transposon_annotation")
                })
    return amr_comp_transposon_annotations

def store_amr_comp_transposon(amr_comp_transposon_annotations, cursor, gbk_files):

    output_paths = [str(f) for f in gbk_files]
    output_path_str = ";".join(output_paths) if gbk_files else ""
    total_updated = 0
    for row in amr_comp_transposon_annotations:
        genome_id = row.get("genome_id")
        gene      = row.get("amr_gene")
        contig_id = row.get("contig_id")
        start     = row.get("start")
        ann       = row.get("composite_transposon_annotation")
        cursor.execute(
            """
            UPDATE annotations
               SET composite_transposon_annotation      = ?,
                   composite_transposon_output_path = ?
             WHERE genome_id = ?
               AND gene_name = ?
               AND contig_ID = ?
               AND start = ?
            """,
            (ann, str(output_path_str), genome_id, gene, contig_id, start)
        )
        total_updated += cursor.rowcount
    return total_updated

def parse_tn3_transposon_annotation(tn3_txt_files):
    """
    Parse Tn3 transposon TXT output files.

    Each file has:
      Line 1: QUERY: <contig_id> <organism> ...
      Then *Candidate N* blocks with a Feature table.
      Feature table columns: Feature, Position, Strand, Length(bp), Type, Positives(%), Coverage(%)

    Extract the family name (Type column) from rows where
    Feature == 'transposase'.

    Returns:
        dict keyed by contig_id, each value is a dict with
        'contig_id' and 'tn3_annotation' (comma-separated string).
    """
    tn3_annotations = {}

    for txt_file in tn3_txt_files:
        try:
            contig_id = None
            transposases = []
            in_feature_table = False

            with open(txt_file) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        in_feature_table = False
                        continue

                    # Extract contig_id from QUERY line
                    if line.startswith("QUERY:"):
                        contig_id = line.split()[1]
                        continue

                    # Skip candidate headers, metadata lines
                    if line.startswith("*Candidate") or line.startswith("Length:") or line.startswith("Distance:") or line.startswith("Order:"):
                        continue

                    # Detect the Feature table header
                    if line.startswith("Feature") and "Position" in line:
                        in_feature_table = True
                        continue

                    # Parse Feature table rows
                    if in_feature_table:
                        parts = line.split()
                        if len(parts) >= 5 and parts[0].lower() == "transposase":
                            # parts: [transposase, start..end, strand, length, type, positives%, coverage%]
                            family = parts[4]
                            transposases.append(family)

                    # Stop parsing at sequence sections
                    if line.startswith(">") or line.startswith("Predicted ORFs"):
                        in_feature_table = False
                        break

            if contig_id and transposases:
                if contig_id not in tn3_annotations:
                    tn3_annotations[contig_id] = {
                        "contig_id": contig_id,
                        "tn3_annotation": []
                    }
                tn3_annotations[contig_id]["tn3_annotation"].extend(transposases)

        except Exception as e:
            print(f"Error parsing {txt_file}: {e}")

    # Deduplicate and convert to comma-separated strings
    for contig_id in tn3_annotations:
        elements = tn3_annotations[contig_id]["tn3_annotation"]
        unique_elements = list(set(elements))
        tn3_annotations[contig_id]["tn3_annotation"] = ", ".join(unique_elements)

    return tn3_annotations

def merge_amr_tn3_transposon(amr_annotations, tn3_annotations):
    """Join AMR with Tn3 by contig_id"""
    amr_tn3_annotations = []
    for annotation in amr_annotations:
        info = tn3_annotations.get(annotation["contig_id"])
        if info:
            amr_tn3_annotations.append({
                "genome_id": annotation["genome_id"],
                "contig_id": annotation["contig_id"],
                "start": annotation["start"],
                "amr_gene": annotation["amr_gene"],
                "tn3_annotation": info.get("tn3_annotation")
            })
    return amr_tn3_annotations

def store_amr_tn3_transposon(amr_tn3_annotations, cursor, gbk_files):
    """Update AMR rows with Tn3 info"""
    output_paths = [str(f) for f in gbk_files]
    output_path_str = ";".join(output_paths) if gbk_files else ""

    total_updated = 0
    for row in amr_tn3_annotations:
        genome_id = row.get("genome_id")
        gene      = row.get("amr_gene")
        contig_id = row.get("contig_id")
        start     = row.get("start")
        ann       = row.get("tn3_annotation")

        cursor.execute(
            """
            UPDATE annotations
               SET tn3_transposon_annotation = ?,
                   tn3_transposon_output_path = ?
             WHERE genome_id = ?
               AND gene_name = ?
               AND contig_id = ?
               AND start = ?
            """,
            (ann, output_path_str, genome_id, gene, contig_id, start)
        )
        total_updated += cursor.rowcount
    return total_updated

def parse_integron_annotation(integron_file_path):
    """
    Parse an IntegronFinder .integrons file.

    Only complete integrons are retained (type == 'complete').

    Returns:
        dict keyed by contig_id (ID_replicon), each value is a list of dicts:
        [{"integron_id": str, "start": int, "end": int}, ...]
    """
    try:
        integron_annotations = {}
        # Temporary structure to collect per-integron boundaries
        # Key: (ID_replicon, ID_integron), Value: {"start": int, "end": int}
        integron_bounds = {}

        with open(integron_file_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                # Skip the header line
                if line.startswith("ID_integron"):
                    continue

                parts = line.split("\t")
                if len(parts) < 11:
                    continue

                integron_type = parts[10].strip()
                if integron_type != "complete":
                    continue

                integron_id = parts[0].strip()
                contig_id = parts[1].strip()
                pos_beg = int(parts[3])
                pos_end = int(parts[4])

                key = (contig_id, integron_id)
                if key not in integron_bounds:
                    integron_bounds[key] = {
                        "start": pos_beg,
                        "end": pos_end
                    }
                else:
                    integron_bounds[key]["start"] = min(integron_bounds[key]["start"], pos_beg)
                    integron_bounds[key]["end"] = max(integron_bounds[key]["end"], pos_end)

        # Convert to the contig_id-keyed structure
        for (contig_id, integron_id), bounds in integron_bounds.items():
            if contig_id not in integron_annotations:
                integron_annotations[contig_id] = []
            integron_annotations[contig_id].append({
                "integron_id": integron_id,
                "start": bounds["start"],
                "end": bounds["end"]
            })

        if not integron_annotations:
            print(f"Note: No complete integrons found in {integron_file_path}")

        return integron_annotations

    except Exception as e:
        print(f"An unexpected error occurred while parsing {integron_file_path}: {e}")
        raise

def merge_amr_integron_annotation(amr_annotations, integron_annotations, 
                                  max_distance=5000, contig_map=None):
    amr_integron_annotations = []
    for annotation in amr_annotations:
        regions = integron_annotations.get(annotation["contig_id"])
        if regions:
            cinfo = (contig_map or {}).get(annotation["contig_id"], {})
            proximal_regions = []
            for r in regions:
                if is_proximal(annotation["start"], annotation["stop"],
                               r["start"], r["end"], max_distance,
                               cinfo.get("length"),
                               cinfo.get("circular", False)):
                    proximal_regions.append(
                        f"{r['integron_id']} ({r['start']}-{r['end']})"
                    )
            if proximal_regions:
                amr_integron_annotations.append({
                    "genome_id": annotation["genome_id"],
                    "contig_id": annotation["contig_id"],
                    "start": annotation["start"],
                    "amr_gene": annotation["amr_gene"],
                    "integron_annotation": ", ".join(proximal_regions)
                })
    return amr_integron_annotations

def store_amr_integron_annotations(amr_integron_annotations, cursor, output_path):
    out_str = str(output_path)
    total_updated = 0
    for row in amr_integron_annotations:
        genome_id = row.get("genome_id")
        gene = row.get("amr_gene")
        contig_id = row.get("contig_id")
        start = row.get("start")
        ann = row.get("integron_annotation")
        cursor.execute(
            """
            UPDATE annotations
               SET integron_annotation = ?,
                   integron_output_path = ?
             WHERE genome_id = ?
               AND gene_name = ?
               AND contig_id = ?
               AND start = ?
            """,
            (ann, out_str, genome_id, gene, contig_id, start)
        )
        total_updated += cursor.rowcount
    return total_updated

def insert_genome(cursor, fasta_name, organism, args_dict, ice_element_metadata=None):
    """
    Insert a single genome and all its annotations into the database.
    Called by main() for standalone use, or by batch_insert.py in a loop.

    Parameters:
        cursor: SQLite cursor (caller owns the connection and commits)
        fasta_name: genome ID string (e.g. 'GCA_037674615.1')
        organism: organism name string (e.g. 'Neisseria_gonorrhoeae')
        args_dict: dict with keys:
            sketch_path, sketch_path_published,
            amrfinder_output, amrfinder_output_published,
            contigs_report_path, contigs_report_path_published,
            filtered_hits_report_path, filtered_hits_report_path_published,
            phage_report_path, phage_report_path_published,
            comp_txt_files, comp_txt_files_published,
            tn3_txt_files, tn3_txt_files_published,
            integron_file, integron_file_published,
            max_distance, gbk_path, iceberg_fasta
        ice_element_metadata: pre-loaded dict from build_ice_element_metadata().
            If None, will be loaded from args_dict["iceberg_fasta"].
    """
    # Create genome entry
    genome_id = create_genome_entry(cursor, fasta_name, organism)

    # Record contig lengths and topology for EVERY genome.
    gbk_path = args_dict.get("gbk_path")
    contig_map = {}
    if gbk_path and Path(gbk_path).exists():
        contig_map = build_contig_map(gbk_path)
        store_contigs(build_contig_map(gbk_path), genome_id, cursor)


    # Store sketch file
    sketch_published = args_dict.get("sketch_path_published")
    if sketch_published:
        store_sketch(sketch_published, cursor)

    # Parse and store AMR annotations
    amr_annotations = []
    amrfinder_output = args_dict.get("amrfinder_output")
    if amrfinder_output and Path(amrfinder_output).exists():
        amr_annotations = parse_amr_annotation(amrfinder_output, genome_id)
        if amr_annotations:
            store_amr_annotation(amr_annotations, cursor, args_dict.get("amrfinder_output_published"))
    else:
        print("Note: No AMR file provided or not found; skipping AMR insert.")

    # Parse and store plasmid annotations
    plasmid_annotations = {}
    contigs_report_path = args_dict.get("contigs_report_path")
    if contigs_report_path and Path(contigs_report_path).exists():
        plasmid_annotations = parse_plasmid_annotation(contigs_report_path)
        if plasmid_annotations:
            plasmid_mge_records = [
                {
                    "genome_id": genome_id,
                    "mge_type": "plasmid",
                    "mge_name": record["plasmid"],
                    "contig_id": record["contig_id"],
                    "start_pos": 1,
                    "end_pos": contig_map.get(record["contig_id"], {}).get("length")
                }
                for record in plasmid_annotations.values()
            ]
            store_mge_elements(plasmid_mge_records, cursor, args_dict.get("contigs_report_path_published"))

    if amr_annotations and plasmid_annotations:
        plasmid_amr_annotation = merge_amr_plasmid_annotation(amr_annotations, plasmid_annotations)
        if plasmid_amr_annotation:
            store_amr_plasmid_annotations(plasmid_amr_annotation, cursor, args_dict.get("contigs_report_path_published"))

    # Parse and store ICE annotations (element-level detection)
    ice_annotations = {}
    filtered_hits_report_path = args_dict.get("filtered_hits_report_path")
    if filtered_hits_report_path and Path(filtered_hits_report_path).exists():
        protein_contig_map = {}
        gbk_path = args_dict.get("gbk_path")
        if gbk_path and Path(gbk_path).exists():
            protein_contig_map = build_protein_contig_map(gbk_path, fasta_name)

        # Use pre-loaded metadata if available, otherwise load it
        if ice_element_metadata is None:
            ice_element_metadata = {}
            iceberg_fasta = args_dict.get("iceberg_fasta")
            if iceberg_fasta and Path(iceberg_fasta).exists():
                ice_element_metadata = build_ice_element_metadata(iceberg_fasta)

        ice_annotations = parse_ice_annotation(
            filtered_hits_report_path, protein_contig_map, ice_element_metadata)

        if ice_annotations:
            ice_mge_records = []
            for contig_id, elements in ice_annotations.items():
                for elem in elements:
                    ice_mge_records.append({
                        "genome_id": genome_id,
                        "mge_type": "ice",
                        "mge_name": elem["ice_label"],
                        "contig_id": contig_id,
                        "start_pos": elem["start"],
                        "end_pos": elem["end"]
                    })
            store_mge_elements(ice_mge_records, cursor, args_dict.get("filtered_hits_report_path_published"))

    max_distance = args_dict.get("max_distance", 5000)

    if amr_annotations and ice_annotations:
        ice_amr_annotations = merge_amr_ice_annotation(
            amr_annotations, ice_annotations, max_distance, contig_map)
        if ice_amr_annotations:
            store_amr_ice_annotations(ice_amr_annotations, cursor, args_dict.get("filtered_hits_report_path_published"))

    # Parse and store prophage annotations
    phage_annotations = {}
    phage_report_path = args_dict.get("phage_report_path")
    if phage_report_path and Path(phage_report_path).exists():
        phage_annotations = parse_phage_annotation(phage_report_path)
        if phage_annotations:
            phage_mge_records = []
            for contig_id, regions in phage_annotations.items():
                for region in regions:
                    phage_mge_records.append({
                        "genome_id": genome_id,
                        "mge_type": "prophage",
                        "mge_name": region["prophage_id"],
                        "contig_id": contig_id,
                        "start_pos": region["start"],
                        "end_pos": region["end"]
                    })
            store_mge_elements(phage_mge_records, cursor, args_dict.get("phage_report_path_published"))

    if amr_annotations and phage_annotations:
        phage_amr_annotation = merge_amr_phage_annotation(amr_annotations, phage_annotations, max_distance, contig_map)
        if phage_amr_annotation:
            store_amr_phage_annotations(phage_amr_annotation, cursor, args_dict.get("phage_report_path_published"))

    # Parse and store composite transposon annotations
    comp_transposon_annotations = {}
    comp_txt_list = []
    comp_txt_files = args_dict.get("comp_txt_files")
    if comp_txt_files:
        for p in comp_txt_files:
            if Path(p).exists():
                comp_txt_list.append(str(p))
    if comp_txt_list:
        comp_transposon_annotations = parse_composite_transposon_annotation(comp_txt_list)
        if comp_transposon_annotations:
            comp_mge_records = [
                {
                    "genome_id": genome_id,
                    "mge_type": "composite_transposon",
                    "mge_name": record["composite_transposon_annotation"],
                    "contig_id": record["contig_id"],
                    "start_pos": record.get("start"),
                    "end_pos": record.get("end")
                }
                for record in comp_transposon_annotations.values()
            ]
            store_mge_elements(comp_mge_records, cursor, args_dict.get("comp_txt_files_published"))

    if amr_annotations and comp_transposon_annotations:
        comp_transposon_amr_annotation = merge_amr_comp_transposon(amr_annotations, comp_transposon_annotations, max_distance, contig_map)
        if comp_transposon_amr_annotation:
            store_amr_comp_transposon(comp_transposon_amr_annotation, cursor, args_dict.get("comp_txt_files_published"))

    # Parse and store Tn3 transposon annotations
    tn3_annotations = {}
    tn3_txt_list = []
    tn3_txt_files = args_dict.get("tn3_txt_files")
    if tn3_txt_files:
        tn3_txt_list = [str(p) for p in tn3_txt_files if Path(p).exists()]
    if tn3_txt_list:
        tn3_annotations = parse_tn3_transposon_annotation(tn3_txt_list)
        if tn3_annotations:
            tn3_mge_records = [
                {
                    "genome_id": genome_id,
                    "mge_type": "tn3_transposon",
                    "mge_name": record["tn3_annotation"],
                    "contig_id": record["contig_id"],
                    "start_pos": None,
                    "end_pos": None
                }
                for record in tn3_annotations.values()
            ]
            store_mge_elements(tn3_mge_records, cursor, args_dict.get("tn3_txt_files_published"))

    if amr_annotations and tn3_annotations:
        tn3_amr_annotation = merge_amr_tn3_transposon(amr_annotations, tn3_annotations)
        if tn3_amr_annotation:
            store_amr_tn3_transposon(tn3_amr_annotation, cursor, args_dict.get("tn3_txt_files_published"))

    # Parse and store integron annotations
    integron_annotations = {}
    integron_file = args_dict.get("integron_file")
    if integron_file and Path(integron_file).exists():
        integron_annotations = parse_integron_annotation(integron_file)
        if integron_annotations:
            integron_mge_records = []
            for contig_id, regions in integron_annotations.items():
                for region in regions:
                    integron_mge_records.append({
                        "genome_id": genome_id,
                        "mge_type": "integron",
                        "mge_name": region["integron_id"],
                        "contig_id": contig_id,
                        "start_pos": region["start"],
                        "end_pos": region["end"]
                    })
            store_mge_elements(integron_mge_records, cursor, args_dict.get("integron_file_published"))

    if amr_annotations and integron_annotations:
        integron_amr_annotation = merge_amr_integron_annotation(amr_annotations, integron_annotations, max_distance, contig_map)
        if integron_amr_annotation:
            store_amr_integron_annotations(integron_amr_annotation, cursor, args_dict.get("integron_file_published"))

def main():
    parser = argparse.ArgumentParser(description='Store sketches and annotations into the ETD DB.')
    parser.add_argument('--db_path', type=Path, default=Path(DATABASE_PATH), help='Path to the SQLite database.')
    parser.add_argument('--fasta_name', type=str, required=True, help='Genome fasta id')
    parser.add_argument('--organism', type=str, default=None, help='Organism name')
    parser.add_argument('--sketch_path', type=Path, help='Path to the all genomes sketch file')
    parser.add_argument('--amrfinder_output', type=Path, help='Path to amrfinderplus TSV file')
    parser.add_argument('--contigs_report_path', type=Path, help='Path to mobsuite contigs_report.txt')
    parser.add_argument('--filtered_hits_report_path', type=Path, default=None, help='Path to to ICE filtered_hits TSV')
    parser.add_argument('--phage_report_path', type=Path, default=None, help='Path to the prophage report TSV')
    parser.add_argument('--comp_txt_files', type=Path, nargs='*', default=None, help='Paths to composite transposon TXT files')
    parser.add_argument('--tn3_txt_files', type=Path, nargs='*', default=None,help='Paths to Tn3 transposon TXT files')
    parser.add_argument('--integron_file', type=Path, default=None, help='Path to IntegronFinder .integrons file')
    parser.add_argument('--max_distance', type=int, default=5000, help='Maximum distance (bp) between AMR gene and MGE element to consider them co-located (default is 5000)')
    parser.add_argument('--gbk_path', type=Path, default=None,
                        help='Path to genome GBK file (for protein-to-contig mapping in ICE analysis)')
    parser.add_argument('--iceberg_fasta', type=Path, default=None,
                        help='Path to ICEberg reference FASTA (for element metadata)')

   
    # Published directory paths (for storing in database)
    parser.add_argument('--sketch_path_published', type=str, help='Published path for sketch file')
    parser.add_argument('--amrfinder_output_published', type=str, help='Published path for AMR output')
    parser.add_argument('--contigs_report_path_published', type=str, help='Published path for contigs report')
    parser.add_argument('--filtered_hits_report_path_published', type=str, default=None, help='Published path for ICE output')
    parser.add_argument('--phage_report_path_published', type=str, default=None, help='Published path for phage output')
    parser.add_argument('--comp_txt_files_published', type=str, nargs='*', default=None, help='Published paths to composite transposon TXT files')
    parser.add_argument('--tn3_txt_files_published', type=str, nargs='*', default=None, help='Published paths to Tn3 transposon TXT files')
    parser.add_argument('--integron_file_published', type=str, default=None, help='Published path for IntegronFinder output')

    args = parser.parse_args()

    # Normalize comp_gbk_files to a list of strings (existing only)
#    comp_txt_list = []
#    if args.comp_txt_files:
#        for p in args.comp_txt_files:
#            #p = Path(p)
#            if p.exists():
#                comp_txt_list.append(str(p))


    db_path = Path(args.db_path)
    init_db(db_path)
    conn = connect(db_path)
    cursor = conn.cursor()

    # Create genome entry
 #   genome_id = create_genome_entry(cursor, args.fasta_name, args.organism)


    # store sketch file
 #   if args.sketch_path and Path(args.sketch_path).exists():
 #       store_sketch(args.sketch_path_published, cursor)

    # parse and store amr_annotations
 #   amr_annotations = []
 #   if args.amrfinder_output and Path(args.amrfinder_output).exists():
 #       amr_annotations = parse_amr_annotation(args.amrfinder_output, genome_id)
 #       if amr_annotations:
 #           store_amr_annotation(amr_annotations, cursor, args.amrfinder_output_published)
 #   else:
 #       print("Note: No AMR file provided or not found; skipping AMR insert.")

    # parse and store plasmid annotations
 #   plasmid_annotations = {}
 #   if args.contigs_report_path and Path(args.contigs_report_path).exists():
 #       plasmid_annotations = parse_plasmid_annotation(args.contigs_report_path)
 #        # Store ALL plasmid elements in mge_elements table
 #       if plasmid_annotations:
 #           plasmid_mge_records = [
 #               {
 #                   "genome_id": genome_id,
 #                   "mge_type": "plasmid",
 #                   "mge_name": record["plasmid"],
 #                   "contig_id": record["contig_id"],
 #                   "start_pos": None,
 #                   "end_pos": None
 #               }
 #               for record in plasmid_annotations.values()
 #           ]
 #           store_mge_elements(plasmid_mge_records, cursor, args.contigs_report_path_published)

    # Merge with AMR only if AMR annotations exist
 #   if amr_annotations and plasmid_annotations:
 #       plasmid_amr_annotation = merge_amr_plasmid_annotation(amr_annotations, plasmid_annotations)
 #       if plasmid_amr_annotation:
 #           store_amr_plasmid_annotations(plasmid_amr_annotation, cursor, args.contigs_report_path_published)

    # parse and store ICE annotations
#    ice_annotations = {}
#    if args.filtered_hits_report_path and Path(args.filtered_hits_report_path).exists():
#        ice_annotations = parse_ice_annotation(args.filtered_hits_report_path)
#        # Store ALL ICE elements in mge_elements table
#        if ice_annotations:
#            ice_mge_records = [
#                {
#                    "genome_id": genome_id,
#                    "mge_type": "ice",
#                    "mge_name": record["ice"],
#                    "contig_id": record["protein_id"],
#                    "start_pos": None,
#                    "end_pos": None
#                }
#                for record in ice_annotations.values()
#            ]
#            store_mge_elements(ice_mge_records, cursor, args.filtered_hits_report_path_published)
#
#    # Merge with AMR only if AMR annotations exist
#    if amr_annotations and ice_annotations:
#        ice_amr_annotations = merge_amr_ice_annotation(amr_annotations, ice_annotations)
#        if ice_amr_annotations:
#            store_amr_ice_annotations(ice_amr_annotations, cursor, args.filtered_hits_report_path_published)

    # parse and store ICE annotations (element-level detection)
#    ice_annotations = {}
#    if args.filtered_hits_report_path and Path(args.filtered_hits_report_path).exists():
        # Build the two lookup maps needed for element-level ICE detection
#        protein_contig_map = {}
#        if args.gbk_path and Path(args.gbk_path).exists():
#            protein_contig_map = build_protein_contig_map(args.gbk_path, args.fasta_name)
#        ice_element_metadata = {}
#        if args.iceberg_fasta and Path(args.iceberg_fasta).exists():
#            ice_element_metadata = build_ice_element_metadata(args.iceberg_fasta)

        # Parse with element-level grouping and 2-category filter
#       ice_annotations = parse_ice_annotation(
#           args.filtered_hits_report_path, protein_contig_map, ice_element_metadata)

        # Store all confirmed ICE elements in mge_elements table
#       if ice_annotations:
#           ice_mge_records = []
#           for contig_id, elements in ice_annotations.items():
#               for elem in elements:
#                   ice_mge_records.append({
#                       "genome_id": genome_id,
#                       "mge_type": "ice",
#                       "mge_name": elem["ice_label"],
#                       "contig_id": contig_id,
#                       "start_pos": elem["start"],
#                       "end_pos": elem["end"]
#                   })
#           store_mge_elements(ice_mge_records, cursor, args.filtered_hits_report_path_published)

    # Merge with AMR only if AMR annotations exist
#   if amr_annotations and ice_annotations:
#      ice_amr_annotations = merge_amr_ice_annotation(
#            amr_annotations, ice_annotations, args.max_distance)
#        if ice_amr_annotations:
#            store_amr_ice_annotations(ice_amr_annotations, cursor, args.filtered_hits_report_path_published)

    # parse and store prophage annotations
#    phage_annotations = {}
#    if args.phage_report_path and Path(args.phage_report_path).exists():
#        phage_annotations = parse_phage_annotation(args.phage_report_path)
        # Store ALL prophage elements in mge_elements table
        # phage_annotations is {contig_id: [{"prophage_id": ..., "start": ..., "end": ...}, ...]}
#        if phage_annotations:
#            phage_mge_records = []
#            for contig_id, regions in phage_annotations.items():
#                for region in regions:
#                    phage_mge_records.append({
#                        "genome_id": genome_id,
#                        "mge_type": "prophage",
#                        "mge_name": region["prophage_id"],
#                        "contig_id": contig_id,
#                        "start_pos": region["start"],
#                        "end_pos": region["end"]
#                    })
#            store_mge_elements(phage_mge_records, cursor, args.phage_report_path_published)

    # Merge with AMR only if AMR annotations exist
#    if amr_annotations and phage_annotations:
#        phage_amr_annotation = merge_amr_phage_annotation(amr_annotations, phage_annotations, args.max_distance)
#        if phage_amr_annotation:
#            store_amr_phage_annotations(phage_amr_annotation, cursor, args.phage_report_path_published)


    # parse and store composite transposon annotations
#    comp_transposon_annotations = {}
#    if comp_txt_list:
#        comp_transposon_annotations = parse_composite_transposon_annotation(comp_txt_list)
        # Store ALL composite transposon elements in mge_elements table
#        if comp_transposon_annotations:
#            comp_mge_records = [
#                {
#                    "genome_id": genome_id,
#                    "mge_type": "composite_transposon",
#                    "mge_name": record["composite_transposon_annotation"],
#                    "contig_id": record["contig_id"],
#                    "start_pos": record.get("start"),
#                    "end_pos": record.get("end")
#                }
#                for record in comp_transposon_annotations.values()
#            ]
#            store_mge_elements(comp_mge_records, cursor, args.comp_txt_files_published)

    # Merge with AMR only if AMR annotations exist
 #   if amr_annotations and comp_transposon_annotations:
 #       comp_transposon_amr_annotation = merge_amr_comp_transposon(amr_annotations, comp_transposon_annotations, args.max_distance)
 #       if comp_transposon_amr_annotation:
 #           store_amr_comp_transposon(comp_transposon_amr_annotation, cursor, args.comp_txt_files_published)

    # parse and store tn3+TA transposon annotations
#    tn3_annotations = {}
#    tn3_txt_list = []
#    if args.tn3_txt_files:
#        tn3_txt_list = [str(p) for p in args.tn3_txt_files if p.exists()]
#    if tn3_txt_list:
#        tn3_annotations = parse_tn3_transposon_annotation(tn3_txt_list)
        # Store ALL Tn3 elements in mge_elements table
#        if tn3_annotations:
#            tn3_mge_records = [
#                {
#                    "genome_id": genome_id,
#                    "mge_type": "tn3_transposon",
#                    "mge_name": record["tn3_annotation"],
#                    "contig_id": record["contig_id"],
#                    "start_pos": None,
#                    "end_pos": None
#                }
#                for record in tn3_annotations.values()
#            ]
#            store_mge_elements(tn3_mge_records, cursor, args.tn3_txt_files_published)

    # Merge with AMR only if AMR annotations exist
#    if amr_annotations and tn3_annotations:
#        tn3_amr_annotation = merge_amr_tn3_transposon(amr_annotations, tn3_annotations)
#        if tn3_amr_annotation:
#            store_amr_tn3_transposon(tn3_amr_annotation, cursor, args.tn3_txt_files_published)


    # parse and store integron annotations
#    integron_annotations = {}
#    if args.integron_file and Path(args.integron_file).exists():
#        integron_annotations = parse_integron_annotation(args.integron_file)
        # Store ALL integron elements in mge_elements table
        # integron_annotations is {contig_id: [{"integron_id": ..., "start": ..., "end": ...}, ...]}
#        if integron_annotations:
#            integron_mge_records = []
#            for contig_id, regions in integron_annotations.items():
#                for region in regions:
#                    integron_mge_records.append({
#                        "genome_id": genome_id,
#                        "mge_type": "integron",
#                        "mge_name": region["integron_id"],
#                        "contig_id": contig_id,
#                        "start_pos": region["start"],
#                        "end_pos": region["end"]
#                    })
#            store_mge_elements(integron_mge_records, cursor, args.integron_file_published)

    # Merge with AMR only if AMR annotations exist
#    if amr_annotations and integron_annotations:
#        integron_amr_annotation = merge_amr_integron_annotation(amr_annotations, integron_annotations, args.max_distance)
#        if integron_amr_annotation:
#            store_amr_integron_annotations(integron_amr_annotation, cursor, args.integron_file_published)

    args_dict = {
        "sketch_path": args.sketch_path,
        "sketch_path_published": args.sketch_path_published,
        "amrfinder_output": args.amrfinder_output,
        "amrfinder_output_published": args.amrfinder_output_published,
        "contigs_report_path": args.contigs_report_path,
        "contigs_report_path_published": args.contigs_report_path_published,
        "filtered_hits_report_path": args.filtered_hits_report_path,
        "filtered_hits_report_path_published": args.filtered_hits_report_path_published,
        "phage_report_path": args.phage_report_path,
        "phage_report_path_published": args.phage_report_path_published,
        "comp_txt_files": args.comp_txt_files,
        "comp_txt_files_published": args.comp_txt_files_published,
        "tn3_txt_files": args.tn3_txt_files,
        "tn3_txt_files_published": args.tn3_txt_files_published,
        "integron_file": args.integron_file,
        "integron_file_published": args.integron_file_published,
        "max_distance": args.max_distance,
        "gbk_path": args.gbk_path,
        "iceberg_fasta": args.iceberg_fasta,
    }

    insert_genome(cursor, args.fasta_name, args.organism, args_dict)

    conn.commit()
    conn.close()


if __name__ == "__main__":
    main()
