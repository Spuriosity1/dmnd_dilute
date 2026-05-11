#!/usr/bin/env python3
import os
import re
import sys
import json
import sqlite3
import argparse
import numpy as np
import io
from multiprocessing import Pool, cpu_count


def adapt_array(arr):
    """Serialize NumPy arrays to BLOBs for SQLite."""
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return sqlite3.Binary(out.read())


def convert_array(text):
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)


sqlite3.register_adapter(np.ndarray, adapt_array)
sqlite3.register_converter("array", convert_array)


# Regex to parse metadata from filename
FILENAME_REGEX = re.compile(
    r'Z1=([\d\-]+,[\d\-]+,[\d\-]+);'
    r'Z2=([\d\-]+,[\d\-]+,[\d\-]+);'
    r'Z3=([\d\-]+,[\d\-]+,[\d\-]+);'
    r'nn(=[\d,]*);'
    r'(p|pZr)=([\d.]+);'
    r'seed=([a-f0-9]+);'
    r'\.stats\.json$'
)


def parse_filename_metadata(filename):
    match = FILENAME_REGEX.match(filename)
    if not match:
        raise ValueError(f"Filename {filename} does not match expected format")

    # Decide which strategy/table based on the prefix
    if match.group(5) == 'p':
        strat = 'random'
        table = 'stats_random'
    else:
        strat = 'Zr'
        table = 'stats_Zr'

    return dict(
        Z1=match.group(1),
        Z2=match.group(2),
        Z3=match.group(3),
        nn=match.group(4)[1:],
        strategy=strat,
        p=float(match.group(6)),
        seed=match.group(7),
        table=table
    )


def parse_stats_file(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)


def create_database(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    for table in ["stats_random", "stats_Zr"]:
        cursor.execute(f'''
            CREATE TABLE IF NOT EXISTS {table} (
                Z1 TEXT,
                Z2 TEXT,
                Z3 TEXT,
                nn TEXT,
                p REAL,
                seed TEXT,
                n_dimers_2 INTEGER,
                n_dimers_4 INTEGER,
                links INTEGER,
                plaqs INTEGER,
                points INTEGER,
                vols INTEGER,
                n_link_parts INTEGER,
                links_wrap BOOLEAN,
                link_cluster_dist BLOB,
                n_plaq_parts INTEGER,
                plaqs_wrap BOOLEAN,
                plaq_cluster_dist BLOB,
                n_vol_parts INTEGER,
                vols_wrap BOOLEAN,
                vol_cluster_dist BLOB
            )
        ''')
    conn.commit()
    return conn


def parse_record(metadata, stats_data):
#    if stats_data is None:
#        raise ValueError("stats_data invalid")
    counts = stats_data.get("counts", {})
    perc = stats_data.get("percolation", {})
    version = stats_data.get("__version__")
    n_dimers = stats_data.get('n_dimers', {'2': 0, '4': 0})
    if n_dimers is None:
        n_dimers = {'2': 0, '4': 0}

    if version is None:
        return None

    return (
        metadata['table'],  # include table name for routing
        metadata['Z1'],
        metadata['Z2'],
        metadata['Z3'],
        metadata['nn'],
        metadata['p'],
        metadata['seed'],
        n_dimers.get('2', 0),
        n_dimers.get('4', 0),
        counts.get('links'),
        counts.get('plaqs'),
        counts.get('points'),
        counts.get('vols'),
        perc.get('n_link_parts'),
        perc.get('links_wrap'),
        np.array(perc.get('link_cluster_dist')),
        perc.get('n_plaq_parts'),
        perc.get('plaqs_wrap'),
        np.array(perc.get('plaq_cluster_dist')),
        perc.get('n_vol_parts'),
        perc.get('vols_wrap'),
        np.array(perc.get('vol_cluster_dist'))
    )


def process_single_file(args):
    """Process a single file - designed to be called in parallel."""
    filename, directory = args
    filepath = os.path.join(directory, filename)
    try:
        metadata = parse_filename_metadata(filename)
        stats_data = parse_stats_file(filepath)
        return parse_record(metadata, stats_data)
    except Exception as e:
        print(f"Error processing {filename}: {e}")
        return None


def insert_chunk(cursor, data_to_insert):
    """Insert records into the proper table based on their source."""
    grouped = {'stats_random': [], 'stats_Zr': []}
    for rec in data_to_insert:
        if rec is None:
            continue
        table = rec[0]
        grouped[table].append(rec[1:])

    for table, records in grouped.items():
        if not records:
            continue
        print(f"\nInserting {len(records)} records into {table}...")
        cursor.executemany(f'''
            INSERT INTO {table} (
                Z1, Z2, Z3, nn, p, seed,
                n_dimers_2, n_dimers_4,
                links, plaqs, points, vols,
                n_link_parts, links_wrap, link_cluster_dist,
                n_plaq_parts, plaqs_wrap, plaq_cluster_dist,
                n_vol_parts, vols_wrap, vol_cluster_dist
            ) VALUES (?, ?, ?, ?, ?, ?,
                      ?, ?, ?, ?, ?, ?,
                      ?, ?, ?,
                      ?, ?, ?,
                      ?, ?, ?)
        ''', records)


def process_directory(directory, db_path, n_workers=None, chunk_size=1000):
    """Process directory with parallel file parsing."""
    if n_workers is None:
        n_workers = cpu_count()

    conn = create_database(db_path)
    cursor = conn.cursor()

    file_list = [f for f in os.listdir(directory) if f.endswith(".stats.json")]
    print(f"Found {len(file_list)} files to process using {n_workers} workers")

    args_list = [(filename, directory) for filename in file_list]

    data_to_insert = []
    with Pool(processes=n_workers) as pool:
        for i, result in enumerate(pool.imap(process_single_file, args_list)):
            if result is not None:
                data_to_insert.append(result)

            if (i + 1) % 100 == 0:
                print(f"Processed {i + 1}/{len(file_list)} ({100.0 * (i + 1) / len(file_list):.1f}%)")

            if len(data_to_insert) >= chunk_size:
                insert_chunk(cursor, data_to_insert)
                conn.commit()
                data_to_insert = []

    if data_to_insert:
        insert_chunk(cursor, data_to_insert)
        conn.commit()

    conn.close()
    print("\nProcessing complete!")
    return file_list


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("DB_REPO", help="Path to the directory containing the json files to combine")
    ap.add_argument("--db", "-o", default="stats.db", help="Output SQLite database file (default: stats.db)")
    ap.add_argument("--cleanup", action='store_true', help="Moves the scanned files to ./trash")

    args = ap.parse_args()
    file_list = process_directory(args.DB_REPO, args.db)

    if args.cleanup:
        print("Remove files? [y/n]")
        choice = input().lower()
        if not choice.startswith('y'):
            sys.exit(0)

        trashdir = os.path.join(args.DB_REPO, 'trash')
        os.makedirs(trashdir, exist_ok=True)

        for f in file_list:
            src = os.path.join(args.DB_REPO, f)
            dst = os.path.join(trashdir, f)
            os.rename(src, dst)

