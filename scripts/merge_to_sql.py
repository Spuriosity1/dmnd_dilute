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
                            r'Z3=([\d\-]+,[\d\-]+,[\d\-]+);' r'nn(=[\d,]*);'
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
        nn= match.group(4)[1:],
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
        link_part_nelem BLOB,
-- plaquette stats
        n_plaq_parts INTEGER,
        plaqs_wrap BOOLEAN,
        plaq_part_nelem BLOB,
-- volume stats
        n_vol_parts INTEGER,
        vols_wrap BOOLEAN,
        vol_part_nelem BLOB
        )
                                   ''')
    conn.commit()
    return conn


from multiprocessing import Pool, cpu_count
import os
import numpy as np

def parse_record(metadata, stats_data):
    counts = stats_data.get("counts", {})
    perc = stats_data.get("percolation", {})
    return (
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
        np.array(perc.get('link_part_nelem')),
        perc.get('n_plaq_parts'),
        perc.get('plaqs_wrap'),
        np.array(perc.get('plaq_part_nelem')),
        perc.get('n_vol_parts'),
        perc.get('vols_wrap'),
        np.array(perc.get('vol_part_nelem'))
    )

def process_single_file(args):
    """Process a single file - designed to be called in parallel"""
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
    print(f"\n\nInserting {len(data_to_insert)} records into SQL...")
    cursor.executemany('''
        INSERT INTO stats (
            Z1, Z2, Z3, nn, p, seed,
            links, plaqs, points, vols,
            n_link_parts, links_wrap, link_part_nelem,
            n_plaq_parts, plaqs_wrap, plaq_part_nelem,
            n_vol_parts, vols_wrap, vol_part_nelem
        ) VALUES (?, ?, ?, ?, ?, ?,
                  ?, ?, ?, ?,
                  ?, ?, ?,
                  ?, ?, ?,
                  ?, ?, ?)
        ''', data_to_insert)

def process_directory(directory, db_path, n_workers=None, chunk_size=1000):
    """
    Process directory with parallel file parsing

    Args:
        directory: Directory containing .stats.json files
        db_path: Path to SQLite database
        n_workers: Number of worker processes (default: cpu_count())
        chunk_size: Number of records to insert at once
    """
    if n_workers is None:
        n_workers = cpu_count()

    conn = create_database(db_path)
    cursor = conn.cursor()

    # Get list of files to process
    file_list = [f for f in os.listdir(directory) if f.endswith(".stats.json")]
    print(f"Found {len(file_list)} files to process using {n_workers} workers")

    # Prepare arguments for parallel processing
    args_list = [(filename, directory) for filename in file_list]

    # Process files in parallel
    data_to_insert = []
    with Pool(processes=n_workers) as pool:
        # Use imap for progress tracking
        for i, result in enumerate(pool.imap(process_single_file, args_list)):
            if result is not None:
                data_to_insert.append(result)

            # Progress indicator
            if (i + 1) % 100 == 0:
                print(f"Processed {i + 1}/{len(file_list)} ({100.0 * (i + 1) / len(file_list):.1f}%)")

            # Insert in chunks
            if len(data_to_insert) >= chunk_size:
                insert_chunk(cursor, data_to_insert)
                data_to_insert = []
                conn.commit()

    # Insert remaining records
    if data_to_insert:
        insert_chunk(cursor, data_to_insert)
        conn.commit()

    conn.close()
    print(f"\nProcessing complete!")



if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "DB_REPO", help="Path to the directory containing the json files to combine")
    ap.add_argument("--db", default="stats.db", help="Output SQLite database file (default: stats.db)")
    ap.add_argument("--cleanup", help="Moves the scanned files to .trash/current date", action='store_true')

    args = ap.parse_args()

    process_directory(args.DB_REPO, args.db)

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

