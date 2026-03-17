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
    parser.add_argument('--excavate_output_dir', type=str, required=True, help='Directory where excavate outputs are located.')
    parser.add_argument('--filtered_vcf_dir', type=str, required=True, help='Directory where guide-filtered vcfs are located.')
    parser.add_argument('--num_strats', type=int, required=True, help='Number of editing strategies to consider for the analysis.')
    parser.add_argument('--strats', nargs="+", default=[], help="Editing strategies that share this analysis pipeline, in list/array form.")
    parser.add_argument('--run_dir', type=str, required=True, help="Higher up run directory, so we can check other editing strategies' completion.")
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
    gene_info=pd.read_csv(args.gene_info, names=['gene','chrom'],sep='\t')
    output_dir=args.output_dir
    excavate_output_dir=args.excavate_output_dir
    filtered_vcf_dir=args.filtered_vcf_dir
    num_strats=args.num_strats
    strats=args.strats
    run_dir=args.run_dir

    # create a dict of genes and their viable genetic variants that can be targeted, based on results from excavate
    gene_targetableSNPs = {}
    gene_numGuides = {}
    for gene in gene_info.gene:
        if os.path.isdir(excavate_output_dir + "/" + gene + '_output'):
            num_files=len(os.listdir(excavate_output_dir + "/" + gene + '_output'))
            if num_files>0:
                cur_guides = pd.read_csv(excavate_output_dir + "/" + gene + '_output/all_guides_summary.csv')
                if sum(cur_guides['no. of guides found (with ref or alt allele)'])>0:
                    gene_targetableSNPs[gene] = list(cur_guides['SNP position'].unique())
                    gene_numGuides[gene] = sum(cur_guides['no. of guides found (with ref or alt allele)'])
    # save targetable snp dictionary for later use
    with open(output_dir + "/metadata/targetable_SNP_locs_by_gene.pkl","wb") as fp:
        pickle.dump(gene_targetableSNPs, fp)

    # vcf formatting info
    sample_list = list(range(1,2549)) # number of people in 1KG
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]
    num_unique_hets_hit=[]
    union_hets_dict={}
    gene_list=[]

    for gene in list(gene_numGuides.keys()):
        gene_list.append(gene)
        # load filtered vcf for the current gene
        vcf = pd.read_table(filtered_vcf_dir + "/" + gene + '_guide_filtered.vcf', comment='#', header=None)
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
    num_df.to_csv(output_dir + '/num_hets_targeted.txt', sep='\t',index=False)
    with open(output_dir + '/unique_hets_hit.pkl', 'wb') as fp:
        pickle.dump(union_hets_dict, fp)

    # # code below is to combine targetable hets results across editing strategies with this analysis pattern
    # counter=0
    # num_files_required=2 # need both the num_hets_targeted.pkl and the unique_hets_hit.pkl to have been created
    # for strat in strats:
    #     if len(os.listdir(run_dir + "/" + strat + "/excavate/het_individuals"))==(num_files_required+1): # +1 because there's also a metadata directory in there
    #         # if the directory is full, it means that the editing strategy has finished running
    #         counter+=1
    # # now check if all editing strategies have completed by checking 'counter'
    # strat_hets={}
    # strat_nums={}
    # if counter==len(strats):
    #     gene_set=set()
    #     # if all have finished running, merge the information into some summary files
    #     for strat in strats:
    #         # load files
    #         with open(run_dir + "/" + strat + "/excavate/het_individuals/unique_hets_hit.pkl", 'rb') as fp:
    #             strat_hets[strat] = pickle.load(fp)
    #         cur_num_df=pd.read_csv(run_dir + "/" + strat + "/excavate/het_individuals/num_hets_targeted.txt", sep='\t')
    #         strat_nums[strat] = cur_num_df
    #         gene_set.update(cur_num_df['gene'])
            
        # # now let's get the new number of unique het individuals captured across the editing strategies
        # gene_list=[]
        # num_unique_hets_all_strats=[]
        # hets_by_gene_across_strats={}
        # for gene in gene_set:
        #     gene_list.append(gene)
        #     cur_hets=set()
        #     # get current gene's hets across strategies
        #     for strat in strats:
        #         if gene in list(strat_hets[strat].keys()):
        #             cur_hets.update(strat_hets[strat][gene])
        #     num_unique_hets_all_strats.append(len(cur_hets))
        #     hets_by_gene_across_strats[gene]=cur_hets
        # # combine into df
        # # super long column name, I know
        # num_df_all_strats=pd.DataFrame({
        #     'gene':gene_list,
        #     'num_unique_hets_captured_across_editing_strategies':num_unique_hets_all_strats
        # })
        # # save
        # num_df_all_strats.to_csv(run_dir + "/summary_files/num_unique_hets_all_strats.csv")
        # with open(run_dir + "/summary_files/unique_hets_across_strats.pkl", 'wb') as fp:
        #     pickle.dump(hets_by_gene_across_strats,fp)




if __name__ == '__main__':
    main()