import os
import re
import json
import sqlite3
import argparse
import json
import sqlite3
import argparse


# Regex to parse metadata from filename
FILENAME_REGEX = re.compile( r'Z1=([\d\-]+,[\d\-]+,[\d\-]+);'
                            r'Z2=([\d\-]+,[\d\-]+,[\d\-]+);'
                            r'Z3=([\d\-]+,[\d\-]+,[\d\-]+);' r'nn=([\d,]+);'
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
        seed TEXT PRIMARY KEY,
        links INTEGER,
        plaqs INTEGER,
        points INTEGER,
        vols INTEGER,
        n_plaq_parts INTEGER,
        n_vol_parts INTEGER
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
            n_plaq_parts, n_vol_parts
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
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
        perc.get('n_plaq_parts'),
        perc.get('n_vol_parts')
    ))

def process_directory(directory, db_path):
    conn = create_database(db_path)
    cursor = conn.cursor()
    file_list = []
    for filename in os.listdir(directory):
        if filename.endswith(".stats.json"):
            filepath = os.path.join(directory, filename)
            file_list.append(filename)
            try:
                metadata = parse_filename_metadata(filename)
                stats_data = parse_stats_file(filepath)
                insert_record(cursor, metadata, stats_data)
            except Exception as e:
                print(f"Error processing {filename}: {e}")
    conn.commit()
    conn.close()

    return file_list




if __name__== "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("DB_REPO", help="Path to the directory containing the json files to combine")
    ap.add_argument("--db", default="stats.db", help="Output SQLite database file (default: stats.db)")
    ap.add_argument("--cleanup", help="Moves the scanned files to .trash/current date", action='store_true')

    args = ap.parse_args()

    file_list= process_directory(args.DB_REPO, args.db)
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

            
        




