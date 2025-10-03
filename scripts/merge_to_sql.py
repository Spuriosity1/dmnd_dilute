import os
import re
import sys
import json
import sqlite3
import argparse
import numpy as np
import io


def adapt_array(arr):
    """
    http://stackoverflow.com/a/31312102/190597 (SoulNibbler)
    """
    out = io.BytesIO()
    np.save(out, arr)
    out.seek(0)
    return sqlite3.Binary(out.read())


def convert_array(text):
    out = io.BytesIO(text)
    out.seek(0)
    return np.load(out)


# Converts np.array to TEXT when inserting
sqlite3.register_adapter(np.ndarray, adapt_array)

# Converts TEXT to np.array when selecting
sqlite3.register_converter("array", convert_array)



# Regex to parse metadata from filename
FILENAME_REGEX = re.compile( r'Z1=([\d\-]+,[\d\-]+,[\d\-]+);'
                            r'Z2=([\d\-]+,[\d\-]+,[\d\-]+);'
                            r'Z3=([\d\-]+,[\d\-]+,[\d\-]+);' r'nn=([\d,]*);'
                            r'p=([\d.]+);' r'seed=([a-f0-9]+);'
                            r'\.stats\.json$')

def parse_filename_metadata(filename):
    match = FILENAME_REGEX.match(filename)
    if not match: 
        raise ValueError(f"Filename {filename} does not match expected format")
    return dict(
        Z1= match.group(1), 
        Z2= match.group(2),
        Z3= match.group(3),
        nn= match.group(4),
        p= float(match.group(5)),
        seed= match.group(6), 
    )

def parse_stats_file(filepath):
    with open(filepath, 'r') as f:
        return json.load(f)


def create_database(db_path):
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS stats (
        Z1 TEXT,
        Z2 TEXT,
        Z3 TEXT,
        nn TEXT,
        p REAL,
        seed TEXT,
        links INTEGER,
        plaqs INTEGER,
        points INTEGER,
        vols INTEGER,
-- link stats
        n_link_parts INTEGER,
        links_wrap BOOLEAN,
        link_cluster_dist BLOB,
-- plaquette stats
        n_plaq_parts INTEGER,
        plaqs_wrap BOOLEAN,
        plaq_cluster_dist BLOB,
-- volume stats
        n_vol_parts INTEGER,
        vols_wrap BOOLEAN,
        vol_cluster_dist BLOB
        )
                                   ''')
    conn.commit()
    return conn

def insert_record(cursor, metadata, stats_data):
    counts = stats_data.get("counts", {})
    perc = stats_data.get("percolation", {})
    cursor.execute('''
        INSERT INTO stats (
            Z1, Z2, Z3, nn, p, seed,
            links, plaqs, points, vols,
            n_link_parts, links_wrap, link_cluster_dist,
            n_plaq_parts, plaqs_wrap, plaq_cluster_dist,
            n_vol_parts, vols_wrap, vol_cluster_dist
        ) VALUES (?, ?, ?, ?, ?, ?,
                  ?, ?, ?, ?,
                  ?, ?, ?,
                  ?, ?, ?,
                  ?, ?, ?)
    ''', (
        metadata['Z1'],
        metadata['Z2'],
        metadata['Z3'],
        metadata['nn'],
        metadata['p'],
        metadata['seed'],
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
    ))

def process_directory(directory, db_path):
    conn = create_database(db_path)
    cursor = conn.cursor()
    file_list = []
    for filename in os.listdir(directory):
        if filename.endswith(".stats.json"):
            filepath = os.path.join(directory, filename)
            try:
                metadata = parse_filename_metadata(filename)
                stats_data = parse_stats_file(filepath)
                insert_record(cursor, metadata, stats_data)
                file_list.append(filename)
            except Exception as e:
                print(f"Error processing {filename}: {e}")
    conn.commit()
    conn.close()

    return file_list


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "DB_REPO", help="Path to the directory containing the json files to combine")
    ap.add_argument("--db", default="stats.db", help="Output SQLite database file (default: stats.db)")
    ap.add_argument("--cleanup", help="Moves the scanned files to .trash/current date", action='store_true')

    args = ap.parse_args()

    file_list = process_directory(args.DB_REPO, args.db)
    print("Processed %d files" % len(file_list))

    if args.cleanup:
        print("Remove files? [y/n]")
        choice = input().lower()
        if not choice.startswith('y'):
            sys.exit(0)

        trashdir=os.path.join(args.DB_REPO, 'trash')
        os.makedirs(trashdir, exist_ok=True)

        for f in file_list:
            filepath = os.path.join(args.DB_REPO, f)
            os.rename(filepath, os.path.join(trashdir, f))

