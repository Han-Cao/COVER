import pytest
import subprocess
import os

import pandas as pd

# Constants
TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(TEST_DIR)
SCRIPT_NAME = 'run_all_het.py'
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
               '-d', os.path.join(TRUTH_DIR, 'Homo_sapiens.GRCh38.110.APP_GFAP.db'),
               '-v', os.path.join(INPUT_DIR, '1kGP_high_coverage_Illumina.GFAP_APP.bcf'),
               '-l', os.path.join(INPUT_DIR, 'APP_GFAP_transcripts.txt'),
               '-p', 'EUR',
               '--top-n-comb', '1000',
               '-o', os.path.join(OUTPUT_DIR, 'APP_GFAP_transcript.EUR')]
    
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

    # read tables
    # we don't compare pair_all since it may have random combinations
    df_test_het_all = pd.read_csv(os.path.join(OUTPUT_DIR, 'APP_GFAP_transcript.EUR.het_freq.all.txt'), sep='\t')
    df_test_het_top = pd.read_csv(os.path.join(OUTPUT_DIR, 'APP_GFAP_transcript.EUR.het_freq.top.txt'), sep='\t')
    #df_test_pair_all = pd.read_csv(os.path.join(OUTPUT_DIR, 'APP_GFAP_transcript.EUR.pair_het_freq.all.txt'), sep='\t')
    df_test_pair_top = pd.read_csv(os.path.join(OUTPUT_DIR, 'APP_GFAP_transcript.EUR.pair_het_freq.top.txt'), sep='\t')

    df_truth_het_all = pd.read_csv(os.path.join(TRUTH_DIR, 'APP_GFAP_transcript.EUR.het_freq.all.txt'), sep='\t')
    df_truth_het_top = pd.read_csv(os.path.join(TRUTH_DIR, 'APP_GFAP_transcript.EUR.het_freq.top.txt'), sep='\t')
    #df_truth_pair_all = pd.read_csv(os.path.join(TRUTH_DIR, 'APP_GFAP_transcript.EUR.pair_het_freq.all.txt'), sep='\t')
    df_truth_pair_top = pd.read_csv(os.path.join(TRUTH_DIR, 'APP_GFAP_transcript.EUR.pair_het_freq.top.txt'), sep='\t')

    # sort table by variants or pairs
    df_test_het_all = df_test_het_all.sort_values(by=['variant1', 'variant2']).reset_index(drop=True)
    df_truth_het_all = df_truth_het_all.sort_values(by=['variant1', 'variant2']).reset_index(drop=True)

    # only compare pair_het_freq for pair het
    df_test_pair_top = df_test_pair_top[['transcript_id', 'gene_id', 'gene_name', 'pair_het_freq']]
    df_truth_pair_top = df_truth_pair_top[['transcript_id', 'gene_id', 'gene_name', 'pair_het_freq']]

    # compare tables
    pd.testing.assert_frame_equal(df_test_het_all, df_truth_het_all)
    pd.testing.assert_frame_equal(df_test_het_top, df_truth_het_top)
    pd.testing.assert_frame_equal(df_test_pair_top, df_truth_pair_top)

    os.remove(os.path.join(OUTPUT_DIR, 'APP_GFAP_transcript.EUR.het_freq.all.txt'))
    os.remove(os.path.join(OUTPUT_DIR, 'APP_GFAP_transcript.EUR.het_freq.top.txt'))
    os.remove(os.path.join(OUTPUT_DIR, 'APP_GFAP_transcript.EUR.pair_het_freq.all.txt'))
    os.remove(os.path.join(OUTPUT_DIR, 'APP_GFAP_transcript.EUR.pair_het_freq.top.txt'))
    if len(os.listdir(OUTPUT_DIR)) == 0:
        os.rmdir(OUTPUT_DIR)