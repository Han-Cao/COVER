#!/usr/bin/env python3

# Parse GTF files and store into a database

import argparse
import logging
import os


import sqlite3
import pandas as pd

# Set up logging
# gtfparse can overwrite logging config, so set it here
logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] - [%(levelname)s]: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

from gtfparse import read_gtf



def parse_gtf(gtf_file: str):
    """Read and parse gtf files"""
    logger = logging.getLogger(__name__)

    # read gtf
    logger.info(f'Parse GTF file: {gtf_file}')
    df_gtf = read_gtf(gtf_file)
    df_gtf = df_gtf.to_pandas()

    # keep chr1-22, X, Y
    df_gtf = df_gtf.loc[df_gtf['seqname'].isin([str(i) for i in range(1, 23)] + ['X', 'Y'])].reset_index(drop=True)

    # transcript table:
    # transcript id, gene id,  gene symbol, coordinate, strand, start codon, stop codon, UTRs, exon list
    df_transcript = df_gtf.loc[df_gtf['feature'] == 'transcript', ['transcript_id', 'gene_id', 'gene_name', 'seqname', 'start', 'end', 'strand']]
    df_transcript = df_transcript.drop_duplicates(subset=['transcript_id']).reset_index(drop=True)
    logger.info(f'Find {df_transcript.shape[0]} unique transcripts')

    # coding region table: transcript id, start codon, stop codon
    df_start_codon = df_gtf.loc[df_gtf['feature'] == 'start_codon', ['transcript_id', 'start', 'end', 'strand']]
    df_end_codon = df_gtf.loc[df_gtf['feature'] == 'stop_codon', ['transcript_id', 'start', 'end', 'strand']]
    # find start and end positions
    # start: 3' end of start codon
    # end: 5' start of stop codon
    df_start_codon['cds_start'] = df_start_codon.apply(lambda x: x['end'] if x['strand'] == '+' else x['start'], axis=1)
    df_end_codon['stop_end'] = df_end_codon.apply(lambda x: x['start'] if x['strand'] == '+' else x['end'], axis=1)

    # if more start codon or end codon is across 2 exons:
    # keep the 3' start codon and 5' stop codon
    # + strand: groupby transcript id and sort, keep the last start codon and first stop codon
    # - strand: groupby transcript id and sort, keep the first start codon and last stop codon
    df_start_codon_plus = df_start_codon.loc[df_start_codon['strand'] == '+'].copy()
    df_start_codon_plus = df_start_codon_plus.sort_values(['transcript_id', 'cds_start']).drop_duplicates(subset=['transcript_id'], keep='last').reset_index(drop=True)
    df_start_codon_minus = df_start_codon.loc[df_start_codon['strand'] == '-'].copy()
    df_start_codon_minus = df_start_codon_minus.sort_values(['transcript_id', 'cds_start']).drop_duplicates(subset=['transcript_id'], keep='first').reset_index(drop=True)
    df_start_codon = pd.concat([df_start_codon_plus, df_start_codon_minus]).reset_index(drop=True)

    df_end_codon_plus = df_end_codon.loc[df_end_codon['strand'] == '+'].copy()
    df_end_codon_plus = df_end_codon_plus.sort_values(['transcript_id', 'stop_end']).drop_duplicates(subset=['transcript_id'], keep='first').reset_index(drop=True)
    df_end_codon_minus = df_end_codon.loc[df_end_codon['strand'] == '-'].copy()
    df_end_codon_minus = df_end_codon_minus.sort_values(['transcript_id', 'stop_end']).drop_duplicates(subset=['transcript_id'], keep='last').reset_index(drop=True)
    df_end_codon = pd.concat([df_end_codon_plus, df_end_codon_minus]).reset_index(drop=True)

    # merge start codon and end codon, only keep transcript with both start and stop codon
    df_cds = df_start_codon[['transcript_id', 'cds_start', 'strand']].merge(df_end_codon[['transcript_id', 'stop_end', 'strand']], how='inner')

    # merge coding region to transcript table
    df_transcript = df_transcript.merge(df_cds, how='inner')
    logger.info(f'Find {df_transcript.shape[0]} transcripts with start codon and stop codon')

    # get exon of each transcript
    df_exon2transcript = df_gtf.loc[df_gtf['feature'] == 'exon', ['transcript_id', 'exon_id']]
    df_exon2transcript = df_exon2transcript.drop_duplicates().reset_index(drop=True)
    # concat exon id to comma separated string for each transcript id
    df_exon2transcript = df_exon2transcript.groupby('transcript_id')['exon_id'].apply(lambda x: ','.join(x)).reset_index()
    df_exon2transcript = df_exon2transcript.rename({'exon_id': 'exon_list'}, axis=1)

    # merge
    df_transcript = df_transcript.merge(df_exon2transcript, how='inner')


    # exon table: 
    # exon id, transcript id, coordinate, strand
    df_exon = df_gtf.loc[df_gtf['feature'] == 'exon', ['exon_id', 'seqname', 'start', 'end', 'strand']]
    df_exon = df_exon.drop_duplicates(subset='exon_id').reset_index(drop=True)

    return df_transcript, df_exon

def create_db(gtf_file, db_file):
    """Create database"""

    logger = logging.getLogger(__name__)
    # Parse GTF file
    df_transcript, df_exon = parse_gtf(gtf_file)

    # Create SQLite3 connection and cursor
    logger.info(f'Creating database: {db_file}')
    if os.path.exists(db_file):
        os.remove(db_file)
    conn = sqlite3.connect(db_file)
    c = conn.cursor()

    # Store df_transcript into database
    c.execute("""
        CREATE TABLE transcripts (
            transcript_id TEXT PRIMARY KEY,
            gene_id TEXT,
            gene_name TEXT,
            seqname TEXT,
            start INTEGER,
            end INTEGER,
            strand TEXT,
            cds_start INTEGER,
            stop_end INTEGER,
            exon_list TEXT
        )""")
    df_transcript.to_sql('transcripts', conn, if_exists='replace', index=False)

    # Create exons_to_transcript table
    c.execute("""
        CREATE TABLE exons (
            exon_id TEXT PRIMARY KEY,
            seqname TEXT,
            start INTEGER,
            end INTEGER,
            strand TEXT
        )""")
    df_exon.to_sql('exons', conn, if_exists='replace', index=False)

    # Commit changes and close connection
    conn.commit()
    conn.close()

def main_gtf2db(args: argparse.Namespace) -> None:
    """Main function"""

    # Create database
    create_db(args.gtf, args.db)

if __name__ == '__main__':
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='Parse GTF files and store into a database')
    parser.add_argument('-g', '--gtf', help='input GTF file')
    parser.add_argument('-d', '--db', help='output SQLite3 database file')
    args = parser.parse_args()

    main_gtf2db(args)