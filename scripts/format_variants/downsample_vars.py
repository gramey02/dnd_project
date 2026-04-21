import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle
from itertools import product,combinations

#parse the arguments passed into the script
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--num_samples', type = int, required = True, help = 'Filepath for the common var dict.')
    parser.add_argument('--output_dir', type=str, required=True, help="output directory for the pipeline run.")
    parser.add_argument('--input_vcf_fp', type=str, required=True, help="Location to save valid snp pairs to.")
    parser.add_argument('--ds_threshold', type=int, required=True, help='number of vars to downsample to')
    args = parser.parse_args()
    return args

def transpose_and_clean(raw_r_tdf):
    r_df = raw_r_tdf.copy() # make a copy
    r_tdf = r_df.transpose() # transpose the copy
    r_tdf.columns = r_tdf.loc["pos", :].values # make the column names the numbers in the position column
    key_df = r_df.loc[
        :, ['chr', 'pos', 'ref', 'alt']
    ].copy()

    clean_r_tdf = r_tdf.iloc[9:, :].copy() # remove some of the initial unneccessary columns

    return clean_r_tdf, key_df
def make_list_of_hets(r_clean_tdf):
    # for each SNP, make a list of heterozygous individuals
    indviduals = r_clean_tdf.index.tolist()
    bool_df = r_clean_tdf.map(lambda x: ((x == "0|1") | (x == "1|0"))).copy()

    snp_het_indviduals_dict = dict() # make a dictionary
    for snp in r_clean_tdf.columns:
        het_individual = bool_df.index[bool_df[snp]].tolist()
        snp_het_indviduals_dict[snp] = het_individual

    return snp_het_indviduals_dict

def assert_unique_and_biallelic_vcf_values(vcf):
    # 1️⃣ remove fully duplicated rows
    vcf = vcf.drop_duplicates()
    
    # 2️⃣ keep only rows where ref and alt are single letters
    vcf = vcf[vcf['ref'].str.len() == 1]
    vcf = vcf[vcf['alt'].str.len() == 1]
    
    # optional: reset index
    vcf = vcf.reset_index(drop=True)
    
    return vcf

def main():
    # parse input arguments
    args = parse_args()
    num_samples=args.num_samples
    output_dir=args.output_dir
    vcf_to_load=args.input_vcf_fp
    ds_vars=args.ds_threshold

    # load some variables useful for vcf loading
    sample_list = list(range(1,num_samples+1)) # number of people in 1KG
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]

    # set genes to downsample
    genes_to_downsample=['DCC']

    # loop through genes
    for gene in genes_to_downsample:
        # load the vcf file for the gene
        vcf = pd.read_table(os.path.join(vcf_to_load, gene + '_CommonVar_filtered.vcf.gz'), comment='#', header=None)
        vcf.columns=cols
        # calculate the most frequently heterozygote frequency for each snp
        vcf = assert_unique_and_biallelic_vcf_values(vcf)
        t_vcf,key_df = transpose_and_clean(vcf)
        total_snps = len(t_vcf.columns)
        het_dict = make_list_of_hets(t_vcf)
        snps=[]
        het_nums=[]
        for snp, hets in het_dict.items():
            snps.append(snp)
            het_nums.append(len(hets))
        paired = sorted(
            zip(snps, het_nums),
            key=lambda x: x[1],
            reverse=True
        )
        # select snps with the top heterozygosity frequency
        if total_snps<ds_vars:
            ds_factor = total_snps
        else:
            ds_factor = ds_vars
        # Take top N specified snp positions
        selected_snps = [pos for pos, count in paired[:ds_factor]]

        # generate a filtered text file so you can filter the vcfs accordingly later
        fp_to_save = os.path.join(output_dir, "excavate/CommonVar_locs", gene + "_CommonVar_locs.txt")
        ds_df = pd.DataFrame({
            'chrom':vcf['chr'].values[0],
            'pos':selected_snps
        })
        ds_df.to_csv(fp_to_save, sep='\t',header=False, index=False)



if __name__=='__main__':
    main()