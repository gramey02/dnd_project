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
DATE = datetime.now().strftime("%Y-%m-%d")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

#parse the arguments passed into the script
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type = str, required = True, help = 'Output directory to save files to.')
    parser.add_argument('--filtered_vcf_dir',type=str,required=True,help='Excavate output directory location.')
    parser.add_argument('--exon_file', type=str,required=True,help="Exon file location.")
    parser.add_argument('--gene', type=str, required=True, help="Curent gene name.")
    parser.add_argument('--num_samples',type=int, required=True, help="Number of samples in population we are determining number of heterozygotes in.")
    parser.add_argument('--excise_entire_gene', type=int, required=True, help="Boolean to tell us if variants should be selected such that the entire gene is excised.")
    parser.add_argument('--valid_pairs_fp', type=str, required=True, help="Path to valid snp pairs that encompass at least one exon for each gene.")
    args = parser.parse_args()
    return args

# load functions for the greedy algorithm
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

def assert_uniq_val_per_row(df, col_name):
    assert (
        df[col_name].nunique() == df[col_name].shape[0]
    ), f"not every row has a unique {col_name} value"

def count_hets_within_region_for_snp_pairs(region_pairs, r_snp_het_indviduals_dict):

    het_counts_in_paired_snps_within_region_df = pd.DataFrame()
    shared_individual_dict = dict()
    for snp1, snp2 in region_pairs:
        n_het_snp1 = len(set(r_snp_het_indviduals_dict[snp1])) # get the number of heterozygous individuals snp 1 has
        n_het_snp2 = len(set(r_snp_het_indviduals_dict[snp2])) # get the number of heterozygous individuals snp 2 has
        #find the overlap between these individuals
        shared_individuals = set(r_snp_het_indviduals_dict[snp1]).intersection(
            set(r_snp_het_indviduals_dict[snp2])
        )
        # make a dictionary where you can lookup the shared individuals for a pair of snps (snp1, snp2)
        shared_individual_dict[(snp1, snp2)] = shared_individuals
        n_shared_individuals = len(shared_individuals)
        n_all_individuals = len(
            set(r_snp_het_indviduals_dict[snp1]).union(
                set(r_snp_het_indviduals_dict[snp2])
            )
        )
        # make a data frame containing the information for a pair of snps (how many each hits individually, and the overlaps/unions of individuals)
        temp_df = pd.DataFrame(
            {
                "snp1": [snp1],
                "snp2": [snp2],
                "n_shared_individuals": [n_shared_individuals],
                "n_het_snp1": [n_het_snp1],
                "n_het_snp2": [n_het_snp2],
                "n_union_individuals": [n_all_individuals],
            }
        )
        #concatenating this newly filled data frame onto the existing one with each iteration
        het_counts_in_paired_snps_within_region_df = (
            pd.concat([het_counts_in_paired_snps_within_region_df, temp_df], ignore_index=True)
        )
        #return the dataframe and the dictionary indicating which snps share which individuals
    return het_counts_in_paired_snps_within_region_df, shared_individual_dict

def encompasses(a,b):
    """Return True if range a is fully encompassed by range b."""
    return (a[0]>=b[0] and a[1]<=b[1])

def pair_encompasses_exon(snp_range,exon_ranges):
    """Return True if snp coordinates fully encompass any of the exons of the gene"""
    encompasses_exon=False
    for exon_range in exon_ranges:
        if encompasses(exon_range,snp_range):
            encompasses_exon=True
            return encompasses_exon
    return encompasses_exon

def is_substring_in_list_loop(substring, string_list):
    for s in string_list:
        if substring in s:
            return True
    return False

def pair_encompasses_exon_ubiquitously(snp_range,exon_info):
    """Return True if snp coordinates fully encompass at least one exon across every transcript of the gene"""
    encompasses_exon=False
    transcripts_where_exon_is_excised=0
    cur_transcript_list=list(exon_info.ensembl_transcript_id.unique())
    total_transcripts=len(cur_transcript_list)

    for transcript in cur_transcript_list:
        transcript_df=exon_info[exon_info.ensembl_transcript_id==transcript]
        cur_exon_ranges=list(set(zip(transcript_df['exon_chrom_start'], transcript_df['exon_chrom_end'])))
        cur_snp_info=pair_encompasses_exon(snp_range,cur_exon_ranges)
        if cur_snp_info==True:
            transcripts_where_exon_is_excised+=1
    # # if the snp range always encompasses an exon, even across all transcripts, return true; else, return False
    # if ((transcripts_where_exon_is_excised<total_transcripts) and (transcripts_where_exon_is_excised>0)):
    #     print(snp_range)
    #     print(transcript)
    return (transcripts_where_exon_is_excised==total_transcripts)

def main():
    # parse input args
    args=parse_args()
    output_dir=args.output_dir
    filtered_vcf_dir=args.filtered_vcf_dir
    exon_file=args.exon_file
    num_samples=args.num_samples
    gene=args.gene
    excise_entire_gene=args.excise_entire_gene
    valid_pairs_fp=args.valid_pairs_fp
    
    # load exon and gene info files
    exon_df=pd.read_csv(exon_file, index_col=0,dtype={'chromosome_name':'str'})

    # get current gene's information
    exon_filt=exon_df[exon_df['hgnc_symbol']==gene][['hgnc_symbol','ensembl_gene_id','chromosome_name','start_position','end_position','strand']].drop_duplicates()
    gene_exons = exon_df[exon_df['hgnc_symbol']==gene]
    gene_exons = gene_exons[gene_exons['transcript_biotype'].isin(['protein_coding','nonsense_mediated_decay', 'non_stop_decay', 'lncRNA', 'miRNA'])]
    gene_start=exon_filt['start_position'].values[0]
    gene_end=exon_filt['end_position'].values[0]

    ##format the vcf column names for ease of use
    sample_list = list(range(1,num_samples+1))
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]

    # create a variable that will store notes on the gene if the greedy algorithm cannot be run for a particular reason
    het_combos_notes=None

    if is_substring_in_list_loop(gene,os.listdir(filtered_vcf_dir)):
        
        # load vcf file for the current gene
        # note that this is the vcf that has already passed allele frequency filtering criteria and that has already been filtered to valid guide positions
        vcf = pd.read_table(filtered_vcf_dir + '/' + gene + '_guide_filtered.vcf', comment='#', header=None)
        vcf.columns = cols

        # create variable that keeps track of whether the criteria to run the greedy algorithm are satisfied
        run_greedy_algorithm=True
            
        # check how many variants lie before and after the gene
        if (excise_entire_gene==1) or (excise_entire_gene==True) or (excise_entire_gene=='True'):
            num_vars_before_start=0
            upstream=[]
            num_vars_after_end=0
            downstream=[]
            for idx,row in vcf.iterrows():
                cur_pos=row.pos
                if cur_pos<gene_start:
                    num_vars_before_start+=1
                    upstream.append(cur_pos)
                if cur_pos>gene_end:
                    num_vars_after_end+=1
                    downstream.append(cur_pos)

            if len(upstream)>0 and len(downstream)>0:
                vcf_filt=vcf[vcf.pos.isin(upstream+downstream)]
            else:
                # code that prevents the running of the rest of the algorithm if there aren't variants that encompass the gene
                run_greedy_algorithm=False
        else:
            vcf_filt=vcf.copy()

        
        if run_greedy_algorithm==True:
            # vcf quality checks
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
            assert_uniq_val_per_row(vcf_filt, "pos") # this function checks that all numbers of the 'pos' column are unique

            # load the gene's valid snp pairs
            with open(os.path.join(valid_pairs_fp, gene+'_valid_snp_pairs.pkl'), 'rb') as fp:
                valid_pairs=pickle.load(fp)

            # Filter this list to just pairs that contain snps that are in the vcf
            vcf_positions=set(vcf_filt['pos'])
            refined_snp_list = {
                pair for pair in valid_pairs
                if pair[0] in vcf_positions and pair[1] in vcf_positions
            }

            # now loop through each passing snp pair and 
            if len(refined_snp_list)>0:
                # proceed with greedy algorithm

                # transpose the vcf data frame and 'clean' it by removing some initial columns that we don't need anymore
                clean_tdf, key_df = transpose_and_clean(vcf_filt)
                
                # output a dictionary where the keys are positions, and the values are lists of samples that are heterozygous for those positions
                snp_het_individuals_dict = make_list_of_hets(clean_tdf)
                
                # get the total number of possible heterozygotes there are
                union_hets = set() # create an empty set
                for snp, het_lists in snp_het_individuals_dict.items():
                    union_hets = set(het_lists).union(union_hets) # combining the snps and samples into a set
                    
                # outputs a dataframe containing information on how many heterozygous people each snp hits indiviudally and how many people pairs of snps share
                (het_counts_snp_pairs_df,shared_ind_dict) = count_hets_within_region_for_snp_pairs(refined_snp_list, snp_het_individuals_dict)
                
                # add a column indicating the proportion of total individuals that the number of shared individuals for each pair of snps represents
                het_counts_snp_pairs_df['prop_shared_individuals'] = het_counts_snp_pairs_df['n_shared_individuals']/num_samples
                het_counts_snp_pairs_df.sort_values("prop_shared_individuals", ascending=False, inplace=True) # sort values in the df putting highest number of individuals hit first
                het_counts_snp_pairs_df.reset_index(drop=True, inplace=True)

                # set up the greedy algorithm below, where we loop through snps, see which one captures the most new heterozygotes, store that one, then remove it from consideration
                top_num_pairs = len(het_counts_snp_pairs_df.snp1.values) # highest number of snp pairs that could possible be investigated
                pair2hets = shared_ind_dict # dictionary: keys = snp pairs, values = sets of heterozygous individuals for the pair
                selected_hets = set()
                selected_pairs = []
                selected_hets_steps = []
                continue_loop=True
                
                # keep the loop running while there are still more snp pairs to hit
                # or if the total number of hets haven't been hit yet
                while (len(selected_pairs)<=top_num_pairs-1) and (len(selected_hets)<len(union_hets)) and (continue_loop==True):
                    best_pair = None
                    best_pair_num_new_hets = 0
                    
                    #one greedy iteration below - for each snp pair and set of hets in the dictionary...
                    for current_pair, current_hets in pair2hets.items():
                        # pick the current snp and all other snps in pair2hets that were previously selected
                        # compare the number of new hets you get to all previously selected hets
                        num_new_hets = len(set(pair2hets[current_pair]) - set(selected_hets))
                        
            
                        # if the number of new hets that this pair adds is greater than previous additions, update "best" parameters
                        if (num_new_hets > best_pair_num_new_hets):
                            best_pair = current_pair
                            best_pair_num_new_hets = num_new_hets
            
                    # update selected variables
                    if best_pair is None:
                        best_pair = current_pair
                    selected_pairs.append(best_pair) # running tally of the best pair selections
                    selected_hets = selected_hets.union(pair2hets[best_pair]) # hets that you have already selected
                    selected_hets_steps.append(len(selected_hets)) # get a stepped readout of how much selected_hets is increasing each time

                    # tell the loop to stop if the same snp pair keeps getting added (i.e., no additional heterozygotes are being added)
                    if len(selected_pairs)>=2:
                        if (selected_pairs[-1]==selected_pairs[-2]):
                            continue_loop=False

                # finally, save the snps that are targeted and the heterozygotes
                with open(output_dir + '/unique_hets_hit_' + gene + '.pkl', 'wb') as fp:
                    pickle.dump(selected_hets,fp)

            else:
                # save information that the greedy algorithm wasn't run because there weren't enough variants that encompass exons
                het_combos_notes="Too few variants encompass exons. Could not run greedy algorithm as a result."
                with open(output_dir + "/" + gene + "no_hets_notes.txt", 'w') as file:
                    file.write(het_combos_notes)
        else:
            # this is executed if the "run_greedy_algorithm" variable is False.
            # Note that the gene didn't have enough variants encompassing it to run the greedy algorithm.
            het_combos_notes="Too few variants upstream or downstream of the gene, and parameters are set so that entire gene is excised."
            with open(output_dir + "/" + gene + "no_hets_notes.txt", 'w') as file:
                file.write(het_combos_notes)
    else:
        # this statement will execute if the gene's file isn't found
        het_combos_notes="Gene's guide-filtered vcf file is not found. Check that vcf filtering occurred properly and didn't exceed run time."
        with open(output_dir + "/" + gene + "no_hets_notes.txt", 'w') as file:
            file.write(het_combos_notes)

if __name__=='__main__':
    main()