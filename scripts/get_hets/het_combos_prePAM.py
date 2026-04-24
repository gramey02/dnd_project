# imports
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle
from datetime import datetime
from itertools import permutations, combinations, product
import itertools

#parse the arguments passed into the script
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type = str, required = True, help = 'Output directory to save files to.')
    parser.add_argument('--filtered_vcf_dir',type=str,required=True,help='Excavate output directory location.')
    parser.add_argument('--valid_pairs_fp', type=str, required=True, help="Valid excision pair file path.")
    parser.add_argument('--num_samples',type=int, required=True, help="Number of samples in population we are determining number of heterozygotes in.")
    parser.add_argument('--gene',type=str, required=True, help="Current gene.")
    args = parser.parse_args()
    return args

# load functions
def make_list_of_hets(r_clean_tdf):
    # for each SNP, make a list of heterozygous individuals
    indviduals = r_clean_tdf.index.tolist()
    bool_df = r_clean_tdf.map(lambda x: ((x == "0|1") | (x == "1|0"))).copy()

    snp_het_indviduals_dict = dict() # make a dictionary
    for snp in r_clean_tdf.columns:
        het_individual = bool_df.index[bool_df[snp]].tolist()
        snp_het_indviduals_dict[snp] = het_individual

    return snp_het_indviduals_dict

def assert_uniq_val_per_row(df, col_name):
    assert (
        df[col_name].nunique() == df[col_name].shape[0]
    ), f"not every row has a unique {col_name} value"

def main():
    # parse input args
    args=parse_args()
    output_dir=args.output_dir
    filtered_vcf_dir=args.filtered_vcf_dir
    num_samples=args.num_samples
    valid_pairs_fp=args.valid_pairs_fp
    gene=args.gene

    ##format the vcf column names for ease of use
    sample_list = list(range(1,num_samples+1))
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]

    het_num_by_gene={}
    hets_by_gene={}

    # load gene's valid excision snp pairs
    with open(valid_pairs_fp + '/' + gene + '_valid_snp_pairs.pkl', 'rb') as fp:
        valid_pairs = pickle.load(fp)

    # load vcf file for the current gene
    # note that this is the vcf that has already passed allele frequency filtering criteria and that has already been filtered to valid guide positions
    vcf = pd.read_table(filtered_vcf_dir + '/' + gene + '_CommonVar_filtered.vcf.gz', compression='gzip', comment='#', header=None)
    vcf.columns = cols

    # vcf quality checks
    vcf_filt=vcf.copy()
    remove_list=[]
    for idx_pos,row_pos in vcf_filt.iterrows():
        # remove indels at the same location as biallelic SNPs that appear as duplicate positions
        if (len(row_pos.alt)>1) or (len(row_pos.ref)>1):
            remove_list.append(1)
        else:
            remove_list.append(0)
        
    vcf_filt=vcf_filt.assign(remove_list=remove_list)
    vcf_filt=vcf_filt[vcf_filt.remove_list==0]
    vcf_filt.drop(labels=['remove_list'],inplace=True,axis=1)
    vcf_filt.drop_duplicates(inplace=True)
    assert_uniq_val_per_row(vcf_filt, "pos") # this function checks that all numbers of the 'pos' column are unique

    # filter the vcf to include only sample information
    sample_cols = vcf_filt.columns[vcf_filt.columns.str.startswith("sample")]
    # boolean DataFrame: True if heterozygous
    is_het = vcf_filt[sample_cols].isin(["0|1", "1|0"])
    # set the snp positions as indices
    is_het_by_pos = is_het.set_axis(vcf_filt["pos"], axis=0)

    remaining_samples = set(sample_cols)

    # run a loop to iteratively determine which people are doubly heterozygous at any of the valid pairs snp locations
    heterozygous_pairs = {}

    for pos1, pos2 in valid_pairs:
        if pos1 not in is_het_by_pos.index or pos2 not in is_het_by_pos.index:
            continue

        s_cols = list(remaining_samples)

        het1 = is_het_by_pos.loc[pos1, s_cols]
        het2 = is_het_by_pos.loc[pos2, s_cols]

        both_het = het1 & het2
        samples = both_het[both_het].index.tolist()

        if samples:
            heterozygous_pairs[(pos1, pos2)] = samples
            remaining_samples -= set(samples)

        if not remaining_samples:
            break

    double_hets=set(sample_cols) - remaining_samples
    total_heterozygous_samples = len(double_hets)
    het_num_by_gene[gene]=total_heterozygous_samples
    hets_by_gene[gene]=double_hets
    
    # save a checkpoint
    with open(os.path.join(output_dir, gene + '_hets_prePAM_checkpoint.pkl'), 'wb') as fp:
        pickle.dump(hets_by_gene, fp)
    with open(os.path.join(output_dir, gene + '_num_hets_by_gene_prePAM_checkpoint.pkl'), 'wb') as fp:
        pickle.dump(het_num_by_gene, fp)


if __name__=='__main__':
    main()