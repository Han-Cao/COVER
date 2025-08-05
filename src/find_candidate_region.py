#!/usr/bin/env python3

# Find candidate deleting region for a given transcript

import argparse
import logging

import pandas as pd
import sqlite3

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] - [%(levelname)s]: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

class Transcript:
    """Class to store transcript information"""
    def __init__(self, info_tuple: tuple) -> None:
        """Initialize from SQL query result"""
        self.id = info_tuple[0]
        self.gene_id = info_tuple[1]
        self.gene_name = info_tuple[2]
        self.seqname = info_tuple[3]
        self.start = info_tuple[4]
        self.end = info_tuple[5]
        self.strand = info_tuple[6]
        self.cds_start = info_tuple[7]
        self.stop_end = info_tuple[8]
        self.exon_list = info_tuple[9].split(",")
        self.exon_count = len(self.exon_list)

        # self.cds_end = self.stop_end - 3 if self.strand == "+" else self.stop_end + 3
        self.seqname = 'chr' + self.seqname if not self.seqname.startswith('chr') else self.seqname

        self.logger = logging.getLogger(__name__)


    def get_exons(self, conn: sqlite3.Connection) -> None:
        """Get exon information from database"""
        exon_query_list = ",".join([f"'{exon_id}'" for exon_id in self.exon_list])
        self.exons = pd.read_sql(f"SELECT * FROM exons WHERE exon_id IN ({exon_query_list})", conn)
        # sort exons from 5' to 3'
        if self.strand == "+":
            self.exons = self.exons.sort_values(by="start").reset_index(drop=True)
        else:
            self.exons = self.exons.sort_values(by="start", ascending=False).reset_index(drop=True)
        
        # get exon length
        self.exons['length'] = self.exons['end'] - self.exons['start'] + 1

        # set exon index as exon no., from 1 to exon_count
        self.exons.index += 1
        
        # get the exon number for start codon and stop codon
        self.start_exon = self.exons.index[(self.exons['start'] <= self.cds_start) & (self.exons['end'] >= self.cds_start)].tolist()
        self.stop_exon = self.exons.index[(self.exons['start'] <= self.stop_end) & (self.exons['end'] >= self.stop_end)].tolist()

        # check if only one exon for start codon and stop codon
        if len(self.start_exon) != 1 or len(self.stop_exon) != 1:
            self.logger.error(f"Transcript {self.id} has more than one exon for start codon or stop codon, please manually select target region")
            exit(1)
        else:
            self.start_exon = self.start_exon[0]
            self.stop_exon = self.stop_exon[0]

        # annotate exons
        self.exons.loc[self.start_exon, 'region'] = 'start codon'
        self.exons.loc[self.stop_exon, 'region'] = 'stop codon'
        if self.start_exon > 1:
            self.exons.loc[self.exons.index < self.start_exon, 'region'] = "5'-UTR"
        if self.stop_exon < self.exon_count:
            self.exons.loc[self.exons.index > self.stop_exon, 'region'] = "3'-UTR"
        if self.stop_exon - self.start_exon > 1:
            self.exons.loc[self.start_exon+1:self.stop_exon-1, 'region'] = "CDS"

        # get bases of frameshift if deletion
        self.exons.loc[self.exons['region'] == 'CDS', 'frameshift'] = self.exons.loc[self.exons['region'] == 'CDS', 'length'] % 3

        # format table
        self.exons['name'] = [f'exon {i}' for i in self.exons.index]
        self.exons = self.exons[['name', 'seqname', 'start', 'end', 'strand', 'length', 'region', 'frameshift', 'exon_id']]
        self.exons['seqname'] = self.exons['seqname'].astype(str).apply(lambda x: 'chr' + x if not x.startswith('chr') else x)

    def get_non_coding(self, max_region: int, splice_donor_len: int, splice_receptor_len: int) -> None:
        """Get non-coding region used for deletion targets"""
        # get coding exons
        df_cds = self.exons.loc[self.exons['region'].isin(['start codon', 'CDS', 'stop codon'])].copy()
        df_cds = df_cds.sort_values(by='start')
        
        # get introns between coding exons
        # intron no.: start codon exon no. ~ stop codon exon no. - 1
        if self.strand == "+":
            intron_no = df_cds.index[:-1]
            df_introns = pd.DataFrame({'seqname': self.seqname,
                                       'start': df_cds['end'].iloc[:-1].values + splice_donor_len + 1,
                                       'end': df_cds['start'].iloc[1:].values - splice_receptor_len - 1,
                                       'strand': self.strand,
                                       'name': [f'intron {i}' for i in intron_no]}, 
                                       index=intron_no)
        else:
            intron_no = df_cds.index[1:]
            df_introns = pd.DataFrame({'seqname': self.seqname,
                                       'start': df_cds['end'].iloc[:-1].values + splice_receptor_len + 1,
                                       'end': df_cds['start'].iloc[1:].values - splice_donor_len - 1,
                                       'strand': self.strand,
                                       'name': [f'intron {i}' for i in intron_no]}, 
                                       index=intron_no)

        # get upstream and downstream region
        # upstream index: 0
        # downstream index: self.exon_count + 1
        # exon index: 1 ~ self.exon_count
        if self.strand == "+":
            df_upstream = pd.DataFrame({'seqname': self.seqname,
                                        'start': self.cds_start - max_region,
                                        'end': self.cds_start - 1,
                                        'strand': self.strand,
                                        'name': "5' region"}, index=[0])
            df_downstream = pd.DataFrame({'seqname': self.seqname,
                                          'start': self.stop_end + 1,
                                          'end': self.stop_end + max_region,
                                          'strand': self.strand,
                                          'name': "3' region"}, index=[self.exon_count + 1])
        else:
            df_upstream = pd.DataFrame({'seqname': self.seqname,
                                        'start': self.cds_start + 1,
                                        'end': self.cds_start + max_region,
                                        'strand': self.strand,
                                        'name': "5' region"}, index=[0])
            df_downstream = pd.DataFrame({'seqname': self.seqname,
                                          'start': self.stop_end - max_region,
                                          'end': self.stop_end - 1,
                                          'strand': self.strand,
                                          'name': "3' region"}, index=[self.exon_count + 1])
        
        # format table
        self.non_coding = pd.concat([df_upstream, df_introns, df_downstream]).sort_index()
        self.non_coding['length'] = self.non_coding['end'] - self.non_coding['start'] + 1
        self.non_coding = self.non_coding[['name', 'seqname', 'start', 'end', 'strand', 'length']]

    def __str__(self) -> str:
        return (f"Transcript ID: {self.id}\n"
                f"Gene ID: {self.gene_id}\n"
                f"Gene Name: {self.gene_name}\n"
                f"Transcript region: {self.seqname}:{self.start}-{self.end} ({self.strand})\n"
                f"Exon count: {self.exon_count}\n"
                f"Start codon exon: {self.start_exon}\n"
                f"Stop codon exon: {self.stop_exon}\n")
        

def read_transcript(conn: sqlite3.Connection) -> list:
    """Get all transcripts in the database"""

    # read all transcript to dataframe
    df_transcript = pd.read_sql('SELECT transcript_id FROM transcripts', conn)

    return df_transcript['transcript_id'].tolist()

def region_dist(x_start: int, x_end: int, y_start: int, y_end: int) -> int:
    """Calculate distance between two genomic regions"""
    if x_start > y_end:
        return x_start - y_end
    elif x_end < y_start:
        return y_start - x_end
    else:
        return 0

    
def find_target_region(x: Transcript, max_deletion: int, n_before_stop: int) -> pd.DataFrame:
    """
    Find target paris of region for deletion
    Creteria:
    1. start codon, or
    2. cause frame shift and exon no. <= stop codon exon no. - N_BEFORE_STOP
    """

    logger = logging.getLogger(__name__)
    lst_target = []
    # find region pairs targeting start codon 
    # upstream region (0 index non-coding region) + any downstream region within max_deletion
    upstream_region = x.non_coding.iloc[0]
    for i in range(1, len(x.non_coding)):
        downstream_region = x.non_coding.iloc[i]
        pair_dist = region_dist(upstream_region['start'], upstream_region['end'], 
                                downstream_region['start'], downstream_region['end'])
        if pair_dist <= max_deletion:
            lst_target.append({'transcript_id': x.id,
                               'gene_id': x.gene_id,
                               'gene_name': x.gene_name,
                               'seqname': x.seqname,
                               'upstream': upstream_region['name'],
                               'upstream_start': upstream_region['start'],
                               'upstream_end': upstream_region['end'],
                               'downstream': downstream_region['name'],
                               'downstream_start': downstream_region['start'],
                               'downstream_end': downstream_region['end'],
                               'distance': pair_dist,
                               'strand': x.strand,
                               'target_exon': f'exon {x.start_exon}',
                               'consequence': 'start loss'})
        else:
            break
    
    # find targetable non-start codon exons
    frameshift_idx = x.exons.index[(x.exons['region'] == 'CDS') & (x.exons.index <= x.stop_exon - n_before_stop) & (x.exons['frameshift'] != 0)]
    if len(lst_target) == 0:
        logger.warning(f"No non-coding regions within {max_deletion}bp can delete start codon")
        logger.info(f"Find {len(frameshift_idx)} target exons to cause frameshift mutation")
    else:
        logger.info(f"Find {len(frameshift_idx) + 1} target exons for knockout (1 start codon, {len(frameshift_idx)} frameshift)")

    # if no targetable exon, return regions targeting start codon
    if len(frameshift_idx) == 0:
        df_target_region = pd.DataFrame(lst_target)
        return df_target_region
    
    # if there are frameshift exons, find region pairs targeting them
    # deletion between intron i and j will delete exon i+1 to exon j
    # upstream: any intron after start codon exon and before stop codon exon - N_BEFORE_STOP
    # downstream: any non-coding region within max_deletion
    # inclusion creteria: sum(frameshift) % 3 != 0
    for upstream_idx in range(x.start_exon, x.stop_exon - n_before_stop): 
        upstream_region = x.non_coding.loc[upstream_idx]
        for downstream_idx in x.non_coding.loc[upstream_idx + 1:].index:
            # if no any frameshift exon in between, skip
            target_exons = x.exons.loc[upstream_idx+1:downstream_idx].index
            if not frameshift_idx.isin(target_exons).any():
                continue
            else:
                # generate target exon string
                if len(target_exons) == 1:
                    target_exons_str = f'exon {target_exons[0]}'
                else:
                    target_exons_str = f'exon {target_exons[0]}-{target_exons[-1]}'
            
            # calculate distance between upstream and downstream regions
            downstream_region = x.non_coding.loc[downstream_idx]
            pair_dist = region_dist(upstream_region['start'], upstream_region['end'], 
                                    downstream_region['start'], downstream_region['end'])
            if pair_dist <= max_deletion:
                # check if sum of frameshift is divisible by 3
                frameshift_sum = x.exons.loc[target_exons, 'frameshift'].sum()
                if frameshift_sum % 3 != 0:
                    lst_target.append({'transcript_id': x.id,
                                       'gene_id': x.gene_id,
                                       'gene_name': x.gene_name,
                                       'seqname': x.seqname,
                                       'upstream': upstream_region['name'],
                                       'upstream_start': upstream_region['start'],
                                       'upstream_end': upstream_region['end'],
                                       'downstream': downstream_region['name'],
                                       'downstream_start': downstream_region['start'],
                                       'downstream_end': downstream_region['end'],
                                       'distance': pair_dist,
                                       'strand': x.strand,
                                       'target_exon': target_exons_str,
                                       'consequence': 'frameshift'})
            # if distance > max_deletion, stop
            else:
                break
    
    df_target_region = pd.DataFrame(lst_target)
    logger.info(f"Find {df_target_region.shape[0]} pairs of candidate regions for allele specific knockout")
    return df_target_region



def query_db(id: str, conn: sqlite3.Connection, max_deletion: int, splice_donor_len: int, splice_receptor_len: int):
    """Query transcript information from database"""
    logger = logging.getLogger(__name__)

    # get cursor
    cur = conn.cursor()

    # Query transcript information
    cur.execute(f"SELECT * FROM transcripts WHERE transcript_id='{id}'")
    res_transcript = cur.fetchone()
    if res_transcript is None:
        logger.error(f"Transcript {id} not found in database, please manually select target region")
        exit(1)
    else:
        transcript = Transcript(res_transcript)
        transcript.get_exons(conn)
        transcript.get_non_coding(max_deletion, splice_donor_len, splice_receptor_len)
    
    cur.close()
    logger.info(f"Retrieve transcript {transcript.id}")
    
    return transcript


def write_output(x: Transcript, df_target_region: pd.DataFrame, out_prefix: str) -> None:
    # write transcript information
    with open(f'{out_prefix}.summary.txt', 'w') as f:
        f.write("# Transcript information:\n")
        f.write(str(x))
        f.write("\n# Exon information:\n")
        f.write(x.exons.to_csv(sep='\t', index=False))
        f.write("\n")
        f.write("\n# Non-coding region information:\n")
        f.write(x.non_coding.to_csv(sep='\t', index=False))
        f.write("\n")

    # write target region information
    df_target_region.to_csv(f'{out_prefix}.candidate_region.txt', sep='\t', index=False)



def main_find_candidate_region(args: argparse.Namespace) -> None:
    """Main function"""

    # Connect to database
    conn = sqlite3.connect(args.db)

    # analysis for list of transcipts
    if args.id_list:
        df_target_region = pd.DataFrame()
        tx_exist = read_transcript(conn)
        tx_missing = []
        with open(args.id_list, 'r') as f:
            for line in f:
                # only analyze transcript in database
                if line.strip() not in tx_exist:
                    tx_missing.append(line.strip())
                else:
                    transcript = query_db(line.strip(), conn, args.max_deletion, args.splice_donor_len, args.splice_receptor_len)
                    df_target_region = pd.concat([df_target_region, find_target_region(transcript, args.max_deletion, args.n_before_stop)])
        # write table
        df_target_region.to_csv(f'{args.output}.candidate_region.txt', sep='\t', index=False)
        # write missing transcript if any
        if len(tx_missing) > 0:
            with open(f'{args.output}.missing_transcript.txt', 'w') as f:
                f.write("\n".join(tx_missing))

    # analysis for single transcript
    elif args.id:
        transcript = query_db(args.id, conn, args.max_deletion, args.splice_donor_len, args.splice_receptor_len)
        df_target_region = find_target_region(transcript, args.max_deletion, args.n_before_stop)
        # write table and summary
        write_output(transcript, df_target_region, args.output)
    else:
        parser.error("Please specify either --id or --id-list")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find candidate region for deletion")
    parser.add_argument("-i", "--id", help="Single transcript ID, will output both candidate region table and summary", type=str)
    parser.add_argument("-l", "--id-list", help="List of transcript IDs to process, will only output candidate region table", type=str)
    parser.add_argument("-d", "--db", help="SQLite3 database file for GTF", required=True, type=str)
    parser.add_argument("-o", "--output", help="Output file prefix", required=True, type=str)

    parser.add_argument("-m", "--max-deletion", help="Maximum deletion size (default: 10000)", type=int, default=10000)
    parser.add_argument("--splice-donor-len", help="Length of splice donor region (default: 10)", type=int, default=10)
    parser.add_argument("--splice-receptor-len", help="Length of splice receptor region (default: 28)", type=int, default=28)
    parser.add_argument("--n-before-stop", help="Minimum number of exons before the stop codon to be considered as target (default: 2)", type=int, default=2)
    args = parser.parse_args()

    main_find_candidate_region(args)