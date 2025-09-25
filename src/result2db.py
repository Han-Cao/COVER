#!/usr/bin/env python3

"""
Script to save analysis results into database
Creates het_freq and pair_het_freq tables and populates them with data from result files
"""

import argparse
import logging
import pandas as pd
import sqlite3
import os

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] - [%(levelname)s]: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

def create_tables(conn: sqlite3.Connection) -> None:
    """Create the het_freq table"""
    cursor = conn.cursor()

    # Create het_freq table
    cursor.execute("""
    CREATE TABLE IF NOT EXISTS het_freq (
        id INTEGER PRIMARY KEY AUTOINCREMENT,
        transcript_id TEXT NOT NULL,
        gene_id TEXT NOT NULL,
        gene_name TEXT,
        variant1 TEXT NOT NULL,
        variant1_region TEXT NOT NULL,
        variant2 TEXT NOT NULL,
        variant2_region TEXT NOT NULL,
        distance INTEGER NOT NULL,
        target TEXT NOT NULL,
        consequence TEXT NOT NULL,
        population TEXT NOT NULL,
        cis_het_freq REAL NOT NULL,
        trans_het_freq REAL NOT NULL,
        max_het_freq REAL NOT NULL,
        target_genotype TEXT NOT NULL,
        created_at TIMESTAMP DEFAULT CURRENT_TIMESTAMP
    )
    """)

    # Create indexes for better query performance
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_het_freq_transcript ON het_freq(transcript_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_het_freq_gene ON het_freq(gene_id)")
    cursor.execute("CREATE INDEX IF NOT EXISTS idx_het_freq_population ON het_freq(population)")

    conn.commit()
    logging.info("Database table created successfully")

def load_het_freq_data(file_path: str) -> pd.DataFrame:
    """Load heterozygous frequency data from file"""
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    df = pd.read_csv(file_path, sep='\t')

    # Validate required columns
    required_columns = [
        'transcript_id', 'gene_id', 'gene_name', 'variant1', 'variant1_region',
        'variant2', 'variant2_region', 'distance', 'target', 'consequence',
        'population', 'cis_het_freq', 'trans_het_freq', 'max_het_freq', 'target_genotype'
    ]

    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing required columns in het_freq file: {missing_columns}")

    logging.info(f"Loaded {len(df)} records from {file_path}")
    return df

def main():
    """Main function"""
    parser = argparse.ArgumentParser(description="Save heterozygous frequency analysis results to database")
    parser.add_argument("-d", "--database", required=True, help="Path to SQLite database file (will be created if it doesn't exist)")
    parser.add_argument("--het-freq-file", required=True, help="Path to heterozygous frequency results file")
    parser.add_argument("--drop-tables", action="store_true", help="Drop existing tables before creating new ones")

    args = parser.parse_args()

    # Check if database file exists before connecting
    db_existed = os.path.exists(args.database)

    # Connect to database (will create if it doesn't exist)
    conn = sqlite3.connect(args.database)

    # Log database creation status
    if db_existed:
        logging.info(f"Using existing database: {args.database}")
    else:
        logging.info(f"Created new database: {args.database}")

    try:
        # Drop tables if requested
        if args.drop_tables:
            cursor = conn.cursor()
            cursor.execute("DROP TABLE IF EXISTS het_freq")
            conn.commit()
            logging.info("Dropped existing tables")

        # Create tables
        create_tables(conn)

        # Process het_freq file
        try:
            het_freq_df = load_het_freq_data(args.het_freq_file)
            het_freq_df.to_sql("het_freq", conn, if_exists="replace", index=False)
        except Exception as e:
            logging.error(f"Error processing het_freq file: {e}")

        logging.info("Database population complete!")

    except Exception as e:
        logging.error(f"Error: {e}")
        conn.rollback()
    finally:
        conn.close()

if __name__ == "__main__":
    main()
