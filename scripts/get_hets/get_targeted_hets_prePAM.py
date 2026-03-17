import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle

#parse the arguments passed into the script
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type = str, required = True, help = 'Output directory for het information.')
    parser.add_argument('--gene_info', type = str, required = True, help = 'Information on genes that have passed filtering criteria so far for the editing strat.')
    parser.add_argument('--filtered_vcf_dir', type=str, required=True, help='Directory where guide-filtered vcfs are located.')
    parser.add_argument('--num_samples', type=int, required=True, help = 'Number of samples in vcf files')
    args = parser.parse_args()
    return args

# functions from Abin that determine num hets
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


def main():
    # parse input args
    args = parse_args()
    gene_info=pd.read_csv(args.gene_info, names=['gene','chrom', 'start_pos', 'end_pos', 'gene_coords', 'io_1', 'io_2'],sep='\t')
    output_dir=args.output_dir
    filtered_vcf_dir=args.filtered_vcf_dir
    num_samples=args.num_samples

    # vcf formatting info
    sample_list = list(range(1,num_samples+1)) # number of people in 1KG
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]
    num_unique_hets_hit=[]
    union_hets_dict={}
    gene_list=[]

    for gene in list(gene_info.gene.unique()):
        gene_list.append(gene)
        # load filtered vcf for the current gene
        vcf = pd.read_table(filtered_vcf_dir + "/" + gene + '_CommonVar_filtered.vcf', comment='#', header=None)
        vcf.columns=cols

        # determine how many people are hit by each snp
        remove=[]
        for idx_pos,row_pos in vcf.iterrows(): # remove duplicate snp positions/indel snp position
            if len(row_pos.alt)>1:
                remove.append(1)
            else:
                remove.append(0)
        vcf['remove'] = remove
        vcf=vcf[vcf['remove']!=1]
        vcf.drop(labels=['remove'],axis=1,inplace=True)
        vcf_clean_tdf, vcf_key_df = transpose_and_clean(vcf)

        # this outputs a dictionary where the keys are positions, and the values are lists of samples that are heterozygous for those positions
        vcf_snp_het_indviduals_dict = make_list_of_hets(vcf_clean_tdf)

        # get a pool of heterozygotes to start with
        union_hets = set() # create an empty set
        for snp, het_lists in vcf_snp_het_indviduals_dict.items():
            union_hets = set(het_lists).union(union_hets) # combining the snps and samples into a set
        num_unique_hets_hit.append(len(union_hets))
        union_hets_dict[gene]=union_hets

    # save the number of individuals the editing strategy in this gene can target
    num_df = pd.DataFrame({
        'gene':gene_list,
        'num_hets_targeted':num_unique_hets_hit
    })
    num_df.to_csv(output_dir + '/num_hets_targeted_prePAM_filter.txt', sep='\t',index=False)
    with open(output_dir + '/unique_hets_hit_prePAM_filter.pkl', 'wb') as fp:
        pickle.dump(union_hets_dict, fp)



if __name__ == '__main__':
    main()