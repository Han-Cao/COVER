#!/usr/bin/env python3

# Calculate co-heterozygous frequency for all transcripts in the database

import argparse
import logging
from multiprocessing import Pool
from functools import partial
import pickle
import os

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] - [%(levelname)s]: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

import pandas as pd
import sqlite3

from find_candidate_region import find_target_region, query_db, Transcript, read_transcript
from calculate_het_freq import calculate_all_het_freq



def transcript_worker(x: Transcript, 
                      max_deletion: int, 
                      n_before_stop: int,
                      vcf_file: str, 
                      pop: str, 
                      maf: float,
                      het: float,
                      pair_per_tx: int,
                      n_pair_max: int,
                      pair_het_cutoff: float,
                      top_n_comb: int) -> tuple:
    """Calculate co-heterozygous frequency for a transcript"""
    # find target region
    df_region = find_target_region(x, max_deletion, n_before_stop)

    # calculate co-heterozygous frequency
    df_het_freq, df_pair_het_freq = calculate_all_het_freq(df_region, vcf_file, pop, maf, het, max_deletion, 
                                                           n_pair_max, pair_het_cutoff, top_n_comb)
    if df_het_freq.shape[0] > pair_per_tx:
        df_het_freq = df_het_freq.sort_values(by='max_het_freq', ascending=False).iloc[:pair_per_tx, :].reset_index(drop=True)

    # add transcript information to the dataframe
    df_het_freq.insert(0, 'transcript_id', x.id)
    df_het_freq.insert(1, 'gene_id', x.gene_id)
    df_het_freq.insert(2, 'gene_name', x.gene_name)

    df_pair_het_freq.insert(0, 'transcript_id', x.id)
    df_pair_het_freq.insert(1, 'gene_id', x.gene_id)
    df_pair_het_freq.insert(2, 'gene_name', x.gene_name)

    return df_het_freq, df_pair_het_freq


def get_max_rows(df: pd.DataFrame, col: str, group: str='transcript_id') -> pd.DataFrame:
    """Get the maximum heterozygous frequency for each transcript"""

    df_max = df.loc[df.groupby(group)[col].idxmax()].reset_index(drop=True)

    return df_max


def main():
    """Main function"""

    logger = logging.getLogger(__name__)

    # parse arguments
    parser = argparse.ArgumentParser(description='Calculate co-heterozygous frequency for all transcripts in the database')
    # required parameters
    parser.add_argument("-d", "--db", help="SQLite3 database file for GTF", required=True)
    parser.add_argument('-v', '--vcf', help='Reference 1000G VCF file', required=True)
    parser.add_argument("-o", "--output", help="Output file prefix", required=True)
    parser.add_argument('-p', '--pop', help='Population', choices=['AFR', 'AMR', 'EAS', 'EUR', 'SAS'], required=True)
    parser.add_argument('--cpu', help='Number of CPUs to use', type=int, default=4)

    # optional parameters for query transcripts
    parser.add_argument('-l', '--id-list', metavar='ID_LIST', help='List of transcript IDs to process')
    parser.add_argument('--save-pickle', help='Save transcript db to pickle file', action='store_true')

    # optional parameters to find candidate region 
    parser.add_argument("--max-deletion", help="Maximum deletion size (default: 10000)", type=int, default=10000)
    parser.add_argument("--splice-donor-len", help="Length of splice donor region (default: 10)", type=int, default=10)
    parser.add_argument("--splice-receptor-len", help="Length of splice receptor region (default: 28)", type=int, default=28)
    parser.add_argument("--n-before-stop", help="Minimum number of exons before the stop codon to be considered as target (default: 2)", type=int, default=2)
    
    # optional parameters to calculate co-heterozygous frequency
    parser.add_argument('--maf', help='MAF cutoff (default: 0.05)', type=float, default=0.05)
    parser.add_argument('--exchet', help='Excess heterozygosity test p-value cutoff (default: 1e-5)', type=float, default=1e-5)
    parser.add_argument('--pair-per-tx', help='Maximum number of variant pairs per transcript to output (default: 1000)', type=int, default=1000)
    parser.add_argument('--n-pair-max', help='Maximum number of variant pairs to be considered for combinations of two variant pairs (default: 200)', type=int, default=200)
    parser.add_argument('--pair-het-cutoff', help='Minimum heterozygous frequency to be considered for combinations of two variant pairs (default: 0.1)', type=float, default=0.1)
    parser.add_argument('--top-n-comb', help='Top N combination of two variant pairs to be output (default: 10)', type=int, default=10)

    args = parser.parse_args()

    # if specified ID list, only work on them
    do_query = True
    if args.id_list:
        with open(args.id_list) as f:
            tx_id_lst = f.read().splitlines()
        conn = sqlite3.connect(args.db)
    else:
        # read all transcripts from database
        db_pickle = f'{args.db}.pickle'
        if os.path.exists(db_pickle):
            logger.info(f'Loading transcript from pickle file: {db_pickle}')
            with open(db_pickle, 'rb') as handle:
                tx_lst = pickle.load(handle)
                do_query = False
        else:
            conn = sqlite3.connect(args.db)
            tx_id_lst = read_transcript(conn)

    # query transcript from database
    if do_query:
        tx_lst = []
        for id in tx_id_lst:
            tx = query_db(id, conn, args.max_deletion, args.splice_donor_len, args.splice_receptor_len)
            tx_lst.append(tx)
        conn.close()
        logger.info(f'Loaded {len(tx_lst)} transcripts from database')
    if args.save_pickle:
        with open(db_pickle, 'wb') as handle:
            pickle.dump(tx_lst, handle)

    # calculate co-heterozygous frequency for each transcript
    # df_het_freq = pd.DataFrame()
    # df_pair_het_freq = pd.DataFrame()

    lst_het_freq = []
    lst_pair_het_freq = []

    # run jobs in parallel
    run_job = partial(transcript_worker, 
                      max_deletion=args.max_deletion, 
                      n_before_stop=args.n_before_stop, 
                      vcf_file=args.vcf,
                      pop=args.pop,
                      maf=args.maf,
                      het=args.exchet,
                      pair_per_tx=args.pair_per_tx,
                      n_pair_max=args.n_pair_max,
                      pair_het_cutoff=args.pair_het_cutoff,
                      top_n_comb=args.top_n_comb)
    
    with Pool(processes=args.cpu) as pool:
        for df_het_freq_tx, df_pair_het_freq_tx in pool.imap_unordered(run_job, tx_lst):
            lst_het_freq.append(df_het_freq_tx)
            lst_pair_het_freq.append(df_pair_het_freq_tx)

    # i = 1
    # for t_id in transcript_list:
    #     logger.info(f'Processing transcript {i}/{n_tx}')
    #     i += 1

    #     # query transcript from database
    #     tx = query_db(t_id, conn, args.max_deletion, args.splice_donor_len, args.splice_receptor_len)

    #     # calculate co-heterozygous frequency
    #     df_het_freq_tx, df_pair_het_freq_tx = transcript_worker(tx, 
    #                                                             args.max_deletion, 
    #                                                             args.n_before_stop, 
    #                                                             args.vcf, 
    #                                                             args.pop, 
    #                                                             args.maf, 
    #                                                             args.pair_per_tx,
    #                                                             args.n_pair_max,
    #                                                             args.pair_het_cutoff,
    #                                                             args.top_n_comb)

    #     # append to dataframe
    #     df_het_freq = pd.concat([df_het_freq, df_het_freq_tx], ignore_index=True)
    #     df_pair_het_freq = pd.concat([df_pair_het_freq, df_pair_het_freq_tx], ignore_index=True)

    # concatenate all dataframes
    df_het_freq = pd.concat(lst_het_freq, ignore_index=True)
    df_pair_het_freq = pd.concat(lst_pair_het_freq, ignore_index=True)
    # get the maximum heterozygous frequency for each transcript
    df_het_freq_max = get_max_rows(df_het_freq, 'max_het_freq')
    df_pair_het_freq_max = get_max_rows(df_pair_het_freq, 'pair_het_freq')

    # write to file
    logger.info(f'Writing to files')
    df_het_freq.to_csv(f'{args.output}.het_freq.all.txt', sep='\t', index=False)
    df_pair_het_freq.to_csv(f'{args.output}.pair_het_freq.all.txt', sep='\t', index=False)
    df_het_freq_max.to_csv(f'{args.output}.het_freq.top.txt', sep='\t', index=False)
    df_pair_het_freq_max.to_csv(f'{args.output}.pair_het_freq.top.txt', sep='\t', index=False)


if __name__ == '__main__':
    main()
