#!/usr/bin/env python3
"""
batch_insert.py – Insert multiple genomes into the ETD database in a single process.

Reads a manifest CSV and calls insert_genome() from aggregate_output.py for each row.
Commits every N genomes (default 500) for efficiency on networked filesystems.
"""

import argparse
import csv
import sys
import time
from pathlib import Path

from aggregate_output import (
    init_db,
    connect,
    insert_genome,
    build_ice_element_metadata,
)


def main():
    parser = argparse.ArgumentParser(
        description="Batch-insert genomes into the ETD database."
    )
    parser.add_argument("--manifest", type=Path, required=True,
                        help="CSV manifest of genomes and annotation paths")
    parser.add_argument("--db_path", type=Path, required=True,
                        help="Path to the SQLite database")
    parser.add_argument("--sketch_path_published", type=str, required=True,
                        help="Published path for the MASH sketch reference")
    parser.add_argument("--iceberg_fasta", type=Path, default=None,
                        help="Path to ICEberg reference FASTA")
    parser.add_argument("--max_distance", type=int, default=5000,
                        help="Max bp distance for AMR-MGE co-location")
    parser.add_argument("--outdir", type=str, required=True,
                        help="Pipeline output directory (for published paths)")
    parser.add_argument("--commit_every", type=int, default=500,
                        help="Commit after every N genomes (default 500)")
    args = parser.parse_args()

    # Initialise database
    init_db(args.db_path)
    conn = connect(args.db_path)
    cursor = conn.cursor()

    # Pre-load ICEberg metadata ONCE
    ice_element_metadata = None
    if args.iceberg_fasta and args.iceberg_fasta.exists():
        print(f"Loading ICEberg metadata from {args.iceberg_fasta} ...")
        ice_element_metadata = build_ice_element_metadata(args.iceberg_fasta)
        print(f"  -> {len(ice_element_metadata)} ICE entries loaded.")

    # Read manifest
    with open(args.manifest, newline="") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    total = len(rows)
    print(f"Manifest contains {total} genomes. Commit interval: {args.commit_every}")

    inserted = 0
    skipped = 0
    n_no_geometry = 0
    t0 = time.time()

    for i, row in enumerate(rows, start=1):
        genome_id = row["genome_id"]
        organism = row.get("organism", "") or None

        # Helper: resolve "NA" / empty -> None, otherwise Path
        def path_or_none(val):
            if not val or val in ("NA", "[]", ""):
                return None
            p = Path(val)
            return p if p.exists() else None

        # Helper: semicolon-delimited multi-file fields -> list of Paths
        def paths_or_none(val):
            if not val or val in ("NA", "[]", ""):
                return None
            parts = [Path(v.strip()) for v in val.split(";") if v.strip()]
            existing = [p for p in parts if p.exists()]
            return existing if existing else None

        # Build the args_dict expected by insert_genome()
        args_dict = {
            "sketch_path_published":              args.sketch_path_published,
            "amrfinder_output":                   path_or_none(row.get("amr_tsv")),
            "amrfinder_output_published":          f"{args.outdir}/amrfinderplus/{genome_id}.tsv",
            "contigs_report_path":                path_or_none(row.get("contigs_report")),
            "contigs_report_path_published":       f"{args.outdir}/mobsuite/{genome_id}_results/contig_report.txt",
            "filtered_hits_report_path":          path_or_none(row.get("ice_hits")),
            "filtered_hits_report_path_published": f"{args.outdir}/filter/{genome_id}_ICEBERG_filtered_hits.tsv",
            "phage_report_path":                  path_or_none(row.get("phage_coords")),
            "phage_report_path_published":         f"{args.outdir}/phispy/{genome_id}_phispy.tsv",
            "comp_txt_files":                     paths_or_none(row.get("comp_txt_files")),
            "comp_txt_files_published":            None,
            "tn3_txt_files":                      paths_or_none(row.get("tn3_txt_files")),
            "tn3_txt_files_published":             None,
            "integron_file":                      path_or_none(row.get("integron_file")),
            "integron_file_published":             f"{args.outdir}/integronfinder/{genome_id}.integrons",
            "geometry_tsv":                       path_or_none(row.get("geometry_tsv")),
            "iceberg_fasta":                      args.iceberg_fasta,
            "max_distance":                       args.max_distance,
        }

        if not args_dict["geometry_tsv"]:
            print(f"WARNING: {genome_id}: geometry TSV not readable "
                  f"({row.get('geometry_tsv')!r}) -- falling back to the GenBank "
                  f"file; if that is also unreadable this genome gets no contigs "
                  f"and no ICE", file=sys.stderr)
            n_no_geometry += 1

        try:
            insert_genome(cursor, genome_id, organism, args_dict,
                          ice_element_metadata=ice_element_metadata)
            inserted += 1
        except Exception as exc:
            import traceback
            traceback.print_exc()
            print(f"ERROR inserting {genome_id}: {exc}", file=sys.stderr)
            skipped += 1
            continue

        # Periodic commit
        if inserted % args.commit_every == 0:
            conn.commit()
            elapsed = time.time() - t0
            rate = inserted / elapsed if elapsed > 0 else 0
            print(f"  Committed {inserted}/{total} genomes  ({rate:.1f} genomes/s)")

    # Final commit
    conn.commit()
    conn.close()

    if n_no_geometry:
        print(f"*** {n_no_geometry}/{total} genomes had no readable geometry TSV. "
              f"Tasks 1-3 and ICE detection are INACTIVE for those genomes. ***",
              file=sys.stderr)
        if n_no_geometry == total:
            sys.exit(1)

    elapsed = time.time() - t0
    print(f"Done. {inserted}/{total} genomes inserted in {elapsed:.1f}s "
          f"({skipped} skipped)")


if __name__ == "__main__":
    main()
