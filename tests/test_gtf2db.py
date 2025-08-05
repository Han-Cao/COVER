import pytest
import subprocess
import os

import pandas as pd
import sqlite3

# Constants
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(TEST_DIR)
SCRIPT_NAME = 'gtf2db.py'
SCRIPT = os.path.join(TEST_DIR, '..', 'src', SCRIPT_NAME)
INPUT_DIR = os.path.join(TEST_DIR, 'input')
TRUTH_DIR = os.path.join(TEST_DIR, 'truth')
OUTPUT_DIR = os.path.join(TEST_DIR, 'output')

# Run command
def run_script() -> None:

    # create output dir if not exist
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

    command = ['python', SCRIPT,
               '-g', os.path.join(INPUT_DIR, 'Homo_sapiens.GRCh38.110.APP_GFAP.gtf'),
               '-d', os.path.join(OUTPUT_DIR, 'Homo_sapiens.GRCh38.110.APP_GFAP.db')]
    
    subprocess.run(command, check=True)


# Test if the command runs successfully
@pytest.mark.order(1)
def test_script_execution() -> None:
    try:
        run_script()  # This will raise CalledProcessError on failure
        assert True
    except subprocess.CalledProcessError as e:
        pytest.fail(f"{SCRIPT_NAME} failed with error: {e}")


# Compare output files with truth
@pytest.mark.order(2)
def test_output() -> None:

    # Extract table from database
    conn_test = sqlite3.connect(os.path.join(OUTPUT_DIR, 'Homo_sapiens.GRCh38.110.APP_GFAP.db'))
    conn_truth = sqlite3.connect(os.path.join(TRUTH_DIR, 'Homo_sapiens.GRCh38.110.APP_GFAP.db'))

    df_test_transcript = pd.read_sql_query("SELECT * FROM transcripts", conn_test)
    df_test_exon = pd.read_sql_query("SELECT * FROM exons", conn_test)
    df_truth_transcript = pd.read_sql_query("SELECT * FROM transcripts", conn_truth)
    df_truth_exon = pd.read_sql_query("SELECT * FROM exons", conn_truth)

    conn_test.close()
    conn_truth.close()

    # Compare tables
    pd.testing.assert_frame_equal(df_test_exon, df_truth_exon, check_dtype=False)
    pd.testing.assert_frame_equal(df_test_transcript, df_truth_transcript, check_dtype=False)

    os.remove(os.path.join(OUTPUT_DIR, 'Homo_sapiens.GRCh38.110.APP_GFAP.db'))
    if len(os.listdir(OUTPUT_DIR)) == 0:
        os.rmdir(OUTPUT_DIR)