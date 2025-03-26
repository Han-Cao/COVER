#!/usr/bin/env python3

# Find deletion target variants within candidate regions and calculate co-heterozygous frequency

import argparse
import logging
import os
from os.path import dirname
from io import StringIO

import pandas as pd
import numpy as np
import pysam
import pysam.bcftools
from itertools import combinations

# Set up logging
logging.basicConfig(level=logging.INFO,
                    format='[%(asctime)s] - [%(levelname)s]: %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

# population sample id files
POP_SAMPLE_FILE = {
    'AFR': os.path.join(dirname(dirname(__file__)), 'data', '1kGP.unrelated.AFR.txt'),
    'AMR': os.path.join(dirname(dirname(__file__)), 'data', '1kGP.unrelated.AMR.txt'),
    'EAS': os.path.join(dirname(dirname(__file__)), 'data', '1kGP.unrelated.EAS.txt'),
    'EUR': os.path.join(dirname(dirname(__file__)), 'data', '1kGP.unrelated.EUR.txt'),
    'SAS': os.path.join(dirname(dirname(__file__)), 'data', '1kGP.unrelated.SAS.txt')
}

# heterozygous genotype code
# homozygous ref/alt: 0
# heterozygous 0|1: 2
# heterozygous 1|0: 3
# combination results: 
# 1. homo + homo = 0 
# 2. homo + hetero <= 3
# 3. cis hetero: 0|1 + 0|1 = 4, 1|0 + 1|0 = 6
# 4. trans hetero: 0|1 + 1|0 = 5
HET_CODE = {
    '0|0': 0,
    '1|1': 0,
    '0|1': 2,
    '1|0': 3
}

HET_CUTOFF = 4
CIS_HET_CODE1 = 4
CIS_HET_CODE2 = 6
TRANS_HET_CODE = 5

class VariantList:
    """Parse VCF genotypes for heterzygous sites"""
    def __init__(self, vcf_file: str, region: str, pop: str, maf: float, het: float) -> None:
        vcf_txt = pysam.bcftools.query("-r", region, 
                                       "-S", POP_SAMPLE_FILE[pop], 
                                       "-i", f"MAF_{pop}_unrel >= {maf} & ExcHet_{pop} >= {het}",
                                       "-f", "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
                                       vcf_file)
        if len(vcf_txt) == 0:
            self.empty = True
            return None
        else:
            self.empty = False
        
        df_vcf = pd.read_csv(StringIO(vcf_txt), sep='\t', header=None)

        # extract variant information
        self.variant = df_vcf.iloc[:, 0] + ':' + df_vcf.iloc[:, 1].astype(str) + ':' + df_vcf.iloc[:, 2] + ':' + df_vcf.iloc[:, 3]
        self.variant = self.variant.to_numpy(dtype=str)
        self.pos = df_vcf.iloc[:, 1].to_numpy(dtype=int)
        self.ref = df_vcf.iloc[:, 2].to_numpy(dtype=str)
        self.alt = df_vcf.iloc[:, 3].to_numpy(dtype=str)

        # extract genotype and conver to het matrix
        self.het_matrix = df_vcf.iloc[:,4:].to_numpy()
        self.het_matrix[self.het_matrix == '0|0'] = HET_CODE['0|0']
        self.het_matrix[self.het_matrix == '1|1'] = HET_CODE['1|1']
        self.het_matrix[self.het_matrix == '0|1'] = HET_CODE['0|1']
        self.het_matrix[self.het_matrix == '1|0'] = HET_CODE['1|0']
        self.het_matrix = self.het_matrix.astype(int)

        self.n_variant = len(self.variant)

    # for loop over variant id and het matrix
    def __iter__(self):
        idx = 0
        while idx < self.n_variant:
            yield self.variant[idx], self.het_matrix[idx]
            idx += 1

    def __len__(self):
        return self.n_variant
    
    def __getitem__(self, idx):
        return self.variant[idx], self.het_matrix[idx]
    
    def get_id(self, idx):
        return self.variant[idx]

    def get_pos(self, idx):
        return self.pos[idx]
    
    def get_ref(self, idx):
        return self.ref[idx]
    
    def get_alt(self, idx):
        return self.alt[idx]
    
    def get_het(self, idx):
        return self.het_matrix[idx]
    
    def get_dist(self, pos: int):
        """Get the distance of all variants to a position"""
        dist = np.abs(self.pos - pos)
        return dist
      
# def var_dist(var1: str, var2: str) -> int:
#     """Calculate distance between two variants with ID format: chr:pos:ref:alt"""
#     var_list1 = var1.split(':')
#     var_list2 = var2.split(':')
#     return abs(int(var_list1[1]) - int(var_list2[1]))


def co_het_freq(x: np.ndarray, het_matrix: np.ndarray) -> tuple:
    """
    Calculate co-heterozygous frequency between SNP x and het matrix (list of SNPs).
    Return cis, trans, and whether most frequent is cis
    """
    # get the number of co-heterozygous samples
    combine_geno = x + het_matrix
    co_het_idx = combine_geno >= HET_CUTOFF
    co_het = co_het_idx.sum(axis=1)
    trans_het_idx = combine_geno == TRANS_HET_CODE
    trans_het = trans_het_idx.sum(axis=1)
    cis_het_idx = co_het_idx != trans_het_idx
    cis_het = co_het - trans_het
    
    # get het frequency
    n_total = x.shape[0]

    cis_het_freq = cis_het / n_total
    trans_het_freq = trans_het / n_total

    # compare cis and trans freq
    is_cis = cis_het >= trans_het

    # get individuals with target genotype (cis_het if is_cis else trans_het)
    # 1 for co-het, 0 for others
    het_ind = get_het_ind(cis_het_idx, trans_het_idx, is_cis)

    return cis_het_freq.round(4), trans_het_freq.round(4), is_cis, het_ind


def get_het_ind(cis_het_idx: np.ndarray, trans_het_idx: np.ndarray, is_cis: np.ndarray) -> np.ndarray:
    """
    Return the het matrix of target het genotype based on is_cis
    cis_het_idx: 2d array of cis het index (N pair * M sample)
    trans_het_idx: 2d array of trans het index (N pair * M sample)
    is_cis: 1d array of whether the most frequent is cis (N pair)
    """
    res_lst = []
    for idx in range(len(is_cis)):
        if is_cis[idx]:
            res_lst.append(cis_het_idx[idx])
        else:
            res_lst.append(trans_het_idx[idx])
    
    return np.array(res_lst)


def target_genotype(ref: str, alt:str, ref_array: np.ndarray, alt_array: np.ndarray, is_cis: np.ndarray) -> np.ndarray:
    """Return the target genotype for one SNP and a list of SNPs"""

    # get alleles for SNP list
    hap1 = np.char.add(ref, np.where(is_cis, ref_array, alt_array))
    hap2 = np.char.add(alt, np.where(is_cis, alt_array, ref_array))

    return np.char.add(np.char.add(hap1, '|'), hap2)


class CoHetMatrix:
    """Store matrix for co-heterozygous sample given a pair of variants"""
    def __init__(self) -> None:
        self.id = []
        self.target = []
        self.cohet_lst = []
    
    def add(self, id: np.ndarray, target: np.ndarray, het_ind: np.ndarray) -> None:
        """
        Add new co-het matrix to the list
        """
        self.id += id.tolist()
        self.target += target.tolist()
        self.cohet_lst += het_ind.tolist() # het_ind: 2d array of het index (N pair * M sample): 1 = het, 0 = hom

    def initialize(self, n_max: int) -> None:
        """
        Convert to matrix for calculating co-het of two pairs of variants
        To calculate pair-het of two pairs of variants, convert the het matrix to hom matrix:
        1 = hom, 0 = het, then np.dot(hom_mat1, hom_mat2.T) denote individuals homozygous for both pairs
        (Note: here "homozygous" means all genotype except the target het genotype, i.e., non cis-het or non trans-het)
        """
        if len(self.id) <= 1:
            self.empty = True
            return None
        else:
            self.empty = False
        
        self.id = np.array(self.id)
        self.target = np.array(self.target)
        self.hom_mat = 1 - np.array(self.cohet_lst)
        self.n_ind = self.hom_mat.shape[1]

        if len(self.id) > n_max:
            # TODO: calculate this twice, can improve
            het_freq = np.sum(self.cohet_lst, axis=1)
            keep_idx = np.argsort(het_freq)[::-1][:n_max]
        
            self.id = self.id[keep_idx]
            self.target = self.target[keep_idx]
            self.hom_mat = self.hom_mat[keep_idx]

    
    def pair_het_table(self) -> pd.DataFrame:
        """Calculate pair het frequency for combinations of two pairs within the matrix"""
        pair_hom = np.dot(self.hom_mat, self.hom_mat.T)
        pair_het_freq = 1 - (pair_hom / self.n_ind) 
        
        # pair_het_freq[i,j] = pair het freq for ids[i] and ids[j]
        # -> upper triangle matrix of pair_het_freq is ids[0,1], ids[0,2] .... ids[0,n] ... ids[n-1, n]
        # -> matched with combinations of ids: combinations(ids, 2)
        comb = list(combinations(self.id, 2))
        pair1 = [x[0] for x in comb]
        pair2 = [x[1] for x in comb]
        het_freq = pair_het_freq[np.triu_indices(len(self.id), k=1)]

        df_pair_het = pd.DataFrame({'pair1': pair1, 'pair2': pair2, 'pair_het_freq': het_freq})

        # map target of pair1 and pair2
        df_pair_target = pd.DataFrame({'pair1': self.id, 'target1': self.target})
        df_pair_het = df_pair_het.merge(df_pair_target, how='left', on='pair1')
        df_pair_target = df_pair_target.rename({'pair1': 'pair2', 'target1': 'target2'}, axis=1)
        df_pair_het = df_pair_het.merge(df_pair_target, how='left', on='pair2')
        df_pair_het = df_pair_het[['pair1', 'target1', 'pair2', 'target2', 'pair_het_freq']]

        return df_pair_het


def pair_het_freq(X: CoHetMatrix, Y: CoHetMatrix) -> np.ndarray:
    """Calculate co-heterozygous frequency for two pairs of variants"""

    pair_hom = np.dot(X.hom_mat, Y.hom_mat.T)
    pair_het_freq = 1 - (pair_hom / X.n_ind) # pair_het_freq[i,j] = pair het freq for X.ids[i] and Y.ids[j]

    return pair_het_freq


def parse_candidate_region(df_region: pd.DataFrame, vcf_file: str, pop: str, maf: float, het: float) -> dict:
    """Parse candidate region dataframe, return a dict of unique VariantList objects"""
    logger = logging.getLogger(__name__)

    res_dict = {}

    for _, row in df_region.iterrows():
        upstream_id = row['upstream']
        
        # add upstream region
        if upstream_id not in res_dict:
            region_str = f"{row['seqname']}:{row['upstream_start']}-{row['upstream_end']}"
            res_dict[upstream_id] = VariantList(vcf_file, region_str, pop, maf, het)
        
        # add downstream region
        downstream_id = row['downstream']
        if downstream_id not in res_dict:
            region_str = f"{row['seqname']}:{row['downstream_start']}-{row['downstream_end']}"
            res_dict[downstream_id] = VariantList(vcf_file, region_str, pop, maf, het)

    logger.info(f'Parse SNP genotypes from {len(res_dict)} unique candidate regions')

    return res_dict


def read_region(region_file: str, idx: list) -> pd.DataFrame:
    """Read candidate region file"""
    logger = logging.getLogger(__name__)
   
    df_region = pd.read_csv(region_file, sep='\t')

    if len(df_region) == 0:
        logger.error(f'No candidate region found in {region_file}')
        exit(1)
    else:
        logger.info(f'Read {df_region.shape[0]} candidate region pairs from {region_file}')

    if len(idx) > 0:
        df_region = df_region.iloc[idx, :].reset_index(drop=True)

    return df_region


def calculate_all_het_freq(df_region: pd.DataFrame, 
                           vcf_file: str, 
                           pop: str, 
                           maf: float,
                           het: float, 
                           max_deletion: int, 
                           n_pair_max: int,
                           pair_het_cutoff: float,
                           top_n_comb: int) -> pd.DataFrame:
    """Calculate co-heterozygous frequency for SNPs in all candidate regions"""
    logger = logging.getLogger(__name__)

    # parse regions
    region_dict = parse_candidate_region(df_region, vcf_file, pop, maf, het)

    # create objects to store results
    df_het_freq = pd.DataFrame()
    co_het_matrix = CoHetMatrix()

    for _, row in df_region.iterrows():
        # make sure region1 is upstream of region2 on plus strand
        if row['upstream_start'] < row['downstream_start']:
            region1 = row['upstream']
            region2 = row['downstream']
        else:
            region1 = row['downstream']
            region2 = row['upstream']
        
        target_id = row['target_exon']
        consequence = row['consequence']

        varlist1 = region_dict[region1]
        varlist2 = region_dict[region2]

        if varlist1.empty or varlist2.empty:
            continue

        # loop over variants in the first region
        for var1_id, var1_het in varlist1:
            # get variants in the second region within deletion distance
            chr, var1_pos, var1_ref, var1_alt = var1_id.split(':')
            var1_pos = int(var1_pos)

            dist = varlist2.get_dist(var1_pos)
            # keep variants within deletion distance
            var2_idx = np.where(dist <= max_deletion)[0]
            if len(var2_idx) == 0:
                continue

            # get co-heterozygous frequency
            cis_freq, trans_freq, is_cis, het_ind = co_het_freq(var1_het, varlist2.get_het(var2_idx))
            max_freq = np.maximum(cis_freq, trans_freq)
            het_geno = target_genotype(var1_ref, var1_alt, varlist2.get_ref(var2_idx), varlist2.get_alt(var2_idx), is_cis)

            # generate dataframe for two SNPs
            df_new = pd.DataFrame({'variant1': var1_id,
                                   'variant1_region': region1,
                                   'variant2': varlist2.get_id(var2_idx),
                                   'variant2_region': region2,
                                   'distance': dist[var2_idx],
                                   'target': target_id,
                                   'consequence': consequence,
                                   'population': pop,
                                   'cis_het_freq': cis_freq,
                                   'trans_het_freq': trans_freq,
                                   'max_het_freq': max_freq,
                                   'target_genotype': het_geno})
            
            # concat to dataframe
            df_het_freq = pd.concat([df_het_freq, df_new], ignore_index=True)

            # add the SNP pair to co-het matrix if max_freq >= pair_het_cutoff
            pair_idx = np.where(max_freq >= pair_het_cutoff)[0]
            if len(pair_idx) == 0:
                continue
            pair_id = np.char.add(f'{var1_id} + ', varlist2.get_id(var2_idx[pair_idx]))
            co_het_matrix.add(pair_id, np.repeat(target_id, pair_id.shape[0]), het_ind[pair_idx])

    # report co-het frequency for one pair of SNPs
    if df_het_freq.empty:
        logger.warning(f'No variant pairs with co-heterozygous frequency >= {pair_het_cutoff}')
        return df_het_freq, pd.DataFrame()

    logger.info(f'Calculate co-heterozygous frequency for {df_het_freq.shape[0]} variant pairs')
    logger.info(f'Maximum co-heterozygous frequency: {df_het_freq["max_het_freq"].max()}')

    # calculate pair-het frequency for two pairs of SNPs
    co_het_matrix.initialize(n_pair_max)
    if co_het_matrix.empty:
        logger.warning(f'Less than 2 variant pairs with co-heterozygous frequency >= {pair_het_cutoff}, skip analysis of SNP pair combinations')
        df_pair_het_freq = pd.DataFrame()
    else:
        logger.info(f'Calculate co-heterozygous frequency for combinations from {len(co_het_matrix.id)} variant pairs')
        df_pair_het_freq = co_het_matrix.pair_het_table()
        df_pair_het_freq = df_pair_het_freq.sort_values('pair_het_freq', ascending=False).reset_index(drop=True).iloc[:top_n_comb, :]

    return df_het_freq, df_pair_het_freq


def main():
    """Main function"""
    parser = argparse.ArgumentParser(description='Calculate co-heterozygous frequency for SNPs in candidate regions')
    parser.add_argument('-r', '--region', help='Candidate region file', required=True)
    parser.add_argument('-v', '--vcf', help='Reference 1000G VCF file', required=True)
    parser.add_argument('-p', '--pop', help='Population', choices=['AFR', 'AMR', 'EAS', 'EUR', 'SAS'], required=True)
    parser.add_argument('-o', '--output', help='Output prefix', required=True)
    parser.add_argument('--index', help='Comma separated index of candidate regions to process (default: all)', type=str, default='all')
    parser.add_argument('-m', '--maf', help='MAF cutoff (default: 0.05)', type=float, default=0.05)
    parser.add_argument('--exchet', help='Excess heterozygosity test p-value cutoff (default: 1e-5)', type=float, default=1e-5)
    parser.add_argument('-d', '--max-deletion', help='Maximum deletion size (default: 10000)', type=int, default=10000)
    parser.add_argument('--n-pair-max', help='Maximum number of variant pairs to be considered for combinations of two variant pairs (default: 200)', type=int, default=200)
    parser.add_argument('--pair-het-cutoff', help='Minimum heterozygous frequency to be considered for combinations of two variant pairs (default: 0.1)', type=float, default=0.1)
    parser.add_argument('--top-n-comb', help='Top N combination of two variant pairs to be output (default: 1000)', type=int, default=1000)

    args = parser.parse_args()

    # parse index string
    if args.index == 'all':
        region_idx = []
    else:
        region_idx = [int(x) for x in args.index.split(',')]

    # read candidate region file
    df_region = read_region(args.region, region_idx)

    # calculate co-heterozygous frequency
    df_het_freq, df_pair_het_freq = calculate_all_het_freq(df_region, args.vcf, args.pop, args.maf, args.exchet,
                                                           args.max_deletion, args.n_pair_max, args.pair_het_cutoff, args.top_n_comb)
    # df_pair_het_freq = df_pair_het_freq.loc[df_pair_het_freq['pair_het_freq'] >= args.pair_het_out_cutoff].reset_index(drop=True)
    # write to file
    df_het_freq.to_csv(f'{args.output}.het_freq.txt', sep='\t', index=False)
    df_pair_het_freq.to_csv(f'{args.output}.pair_het_freq.txt', sep='\t', index=False)

if __name__ == '__main__':
    main()