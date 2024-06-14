"""
Title: AMotifFinder
Description: A program to analyze amino acid sequences that are read from a FASTA file. 
            The program identifies the occurring sub-sequences (sequence motifs) in the file of a certain length, 
            and provides the number and location of them in the protein sequences.
Author: Ala Mańkowska
Email: avrmankowska@gmail.com
Date: 2024-06-14
Version: 1.0
"""
import sys
import uuid
from datetime import datetime
import atexit
from Bio import SeqIO
import csv
import os
from collections import defaultdict
import psutil
import logging
import sqlite3

amino_acids = 'ACDEFGHIKLMNPQRSTVWYBXZ*U'   #IUPAC
aa_to_bits = {aa: idx for idx, aa in enumerate(amino_acids)}
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

db_name = f'{uuid.uuid4()}.db'

def remove_db_file(db_path):
    if os.path.exists(db_path):
        os.remove(db_path)
        logging.info(f"Database file {db_path} removed.")

def encode_amino_acid_sequence(sequence):
    encoded = 0
    for amino_acid in sequence:
        encoded = (encoded << 5) | aa_to_bits[amino_acid]
    return encoded


def decode_amino_acid_sequence(encoded, length):
    decoded = ''
    mask = (1 << 5) - 1 #to extract 5 bits
    for _ in range(length):
        aa_idx = encoded & mask

        decoded = amino_acids[aa_idx] + decoded
        encoded >>= 5
    return decoded


def init_db(subfolder_path):
    conn = sqlite3.connect(os.path.join(subfolder_path, db_name))
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS subsequences (
            length INTEGER,
            subsequence TEXT,
            count INTEGER,
            sequences TEXT
        )
    ''')
    conn.commit()
    conn.execute('PRAGMA journal_mode=WAL;')
    conn.execute('PRAGMA synchronous=NORMAL;')
    return conn


def save_to_db(conn, aggregated_data):
    logging.info(f"Saving {len(aggregated_data)} subsequences to database.")
    c = conn.cursor()
    prepared_data = [
        (length, subseq, info['Count'], ';'.join(info['Sequences']))
        for (length, subseq), info in aggregated_data.items()
    ]
    c.executemany('''
        INSERT INTO subsequences (length, subsequence, count, sequences)
        VALUES (?, ?, ?, ?)
    ''', prepared_data)
    conn.commit()


def count_subsequences(args):
    sequences, min_len, max_len, folder_name = args
    logging.info(f"Processing subsequences from {min_len} to {max_len} in provided sequences.")
    conn = init_db(folder_name)
    subseq_data = defaultdict(lambda: defaultdict(lambda: {'Count': 0, 'Sequences': set()}))
    batch_size = 10000  #Adjust batch size based on your memory capacity and dataset size

    for seq_record in sequences:
        seq_id = seq_record.id.split('|')[1]
        sequence = str(seq_record.seq)
        checked = set()

        for length in range(min_len, max_len + 1):
            for start in range(len(sequence) - length + 1):
                subseq = sequence[start:start + length]

                if subseq in checked:
                    continue

                count = sequence.count(subseq)
                checked.add(subseq)
                subseq_data[length][subseq]['Count'] += count
                subseq_data[length][subseq]['Sequences'].add(seq_id)

                if sum(len(subseqs) for subseqs in subseq_data.values()) >= batch_size:
                    save_to_db(conn,
                               {(length, subseq): data for length, subseqs in subseq_data.items() for subseq, data in
                                subseqs.items()})
                    subseq_data.clear()

    if any(subseq_data.values()):
        save_to_db(conn, {(length, subseq): data for length, subseqs in subseq_data.items() for subseq, data in
                          subseqs.items()})

    conn.close()
    logging.info(f"Finished processing for sequences from {min_len} to {max_len}.")


def export_data_to_csv(length, output_folder):
    logging.info(f"Exporting data to CSV files in folder: {output_folder}")
    if not os.path.exists(output_folder):
        logging.info(f"Creating directory {output_folder}.")
        os.makedirs(output_folder)

    conn = init_db(output_folder)
    cursor = conn.cursor()
    query = '''
        SELECT subsequence, SUM(count) AS TotalCount, GROUP_CONCAT(sequences, ';') AS AllSequences
        FROM subsequences
        WHERE length = ?
        GROUP BY subsequence
    '''
    cursor.execute(query, (length,))
    with open(os.path.join(output_folder, f"{length}.csv"), 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Subsequence', 'Count', 'Sequences'])
        for row in cursor:
            writer.writerow(row)
    cursor.close()
    conn.close()
    logging.info(f"Data exported to CSV files in folder: {output_folder}.")


def print_memory_usage():
    process = psutil.Process(os.getpid())
    logging.info(f"Memory used: {process.memory_info().rss / 1024 ** 2} MB")


if __name__ == "__main__":
    logging.info("Starting the sequence analysis script.")
    if len(sys.argv) < 4:
        print("Usage: python script.py <file_path> <min_len> <max_len>")
        logging.error("Incorrect number of arguments. Exiting.")
        sys.exit(1)

    file_path = sys.argv[1]
    min_len = int(sys.argv[2])
    max_len = int(sys.argv[3])

    try:
        sequences = list(SeqIO.parse(file_path, "fasta"))
        logging.info(f"Loaded {len(sequences)} sequences from {file_path}.")
    except FileNotFoundError:
        print("The FASTA file could not be read.")
        logging.error(f"File {file_path} not found. Exiting.")
        sys.exit(1)

    folder_name = "subsequences_results"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        logging.info(f"Created directory {folder_name}.")

    atexit.register(remove_db_file, os.path.join(folder_name, db_name))

    start = datetime.now()
    count_subsequences((sequences, min_len, max_len, folder_name))
    for length in range(min_len, max_len + 1):
        export_data_to_csv(length, folder_name)

    logging.info(f"Script completed. Total time: {datetime.now() - start}.")
    print_memory_usage()
