# imports
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle

#parse the arguments passed into the script
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type = str, required = True, help = 'Output directory to save files to.')
    parser.add_argument('--exon_file',type=str,required=True,help="Filename for list of chromosomes and genes of interest.")
    parser.add_argument('--af_limit',type=float,required=True,help="Allele frequency threshold.")
    parser.add_argument('--gene_info', type=str,required=True, help="file to the excision windows for each gene.")
    parser.add_argument('--nearby_gene_filter', type=str,required=True, help="Boolean to tell us if excision windows were filtered based on other genes in the window.")
    args = parser.parse_args()
    return args

def encompasses(a,b):
    """Return True if range a is fully encompassed by range b."""
    return (a[0]>=b[0] and a[1]<=b[1])

def main():
    # parse args
    args=parse_args()
    output_dir=args.output_dir
    exon_file=args.exon_file
    af_limit=args.af_limit
    gene_df=pd.read_csv(args.gene_info, sep='\t')
    nearby_gene_filter=args.nearby_gene_filter

    # load exon file
    exon_df=pd.read_csv(exon_file,index_col=0,dtype={'chromosome_name':'str'})
    exon_filt=exon_df[['hgnc_symbol','chromosome_name','start_position','end_position']].drop_duplicates()

    # create vcf dict
    vcf_dict={}
    chroms = gene_df.chrom.unique()
    for chrom in chroms:
        af_filename = '/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/TGP_chr' + str(chrom) + '_afs.txt'
        cur_chrom_TGP_afs = pd.read_csv(af_filename, sep=' ', names = ['chrom', 'pos', 'ref', 'alt', 'ac', 'an', 'af', 'afr_af', 'amr_af', 'eas_af', 'eur_af', 'sas_af'])
        cur_chrom_TGP_afs=cur_chrom_TGP_afs[(cur_chrom_TGP_afs.af>=af_limit) & (cur_chrom_TGP_afs.af<=1-af_limit)]
        # get a list of common variant positions in and around the gene
        vcf_dict[chrom] = cur_chrom_TGP_afs[['chrom','pos','ref','alt', 'af']]

    # now for each gene's editing window, find the variants in that
    common_var_info={}
    num_common_vars_in_window=[]
    gene_list=[]
    for idx,row in gene_df.iterrows():
        cur_gene=row.hgnc_symbol
        cur_chrom=row.chrom
        gene_list.append(cur_gene)
        cur_gene_common_var_info=[]
        if (row.good_excision_candidate==True) or (row.good_excision_candidate==1):
            cur_excision_start=row.new_excision_start
            cur_excision_end=row.new_excision_end

            # get the vcf file specific to this gene's window
            cur_vcf=vcf_dict[cur_chrom]
            cur_vcf_filt = cur_vcf[cur_vcf.pos>=cur_excision_start]
            cur_vcf_filt = cur_vcf_filt[cur_vcf_filt.pos<=cur_excision_end]
            for idx_pos,row_pos in cur_vcf_filt.iterrows():
                cur_gene_common_var_info.append([row_pos.pos,row_pos.af])
            
            num_common_vars_in_window.append(len(cur_gene_common_var_info))
            common_var_info[cur_gene]=cur_gene_common_var_info
        else:
            common_var_info[cur_gene]=cur_gene_common_var_info
            num_common_vars_in_window.append(None)

    # combine information into data frame
    excision_cv_df=pd.DataFrame({
        'hgnc_symbol':gene_list,
        'num_common_vars_in_excision_window':num_common_vars_in_window
    })
    if 'good_excision_candidate' in list(gene_df.columns):
        excision_cv_df['good_excision_candidate']=list(gene_df.good_excision_candidate)
    excision_df_merged = excision_cv_df.merge(exon_filt, on='hgnc_symbol',how='left')
    
    # check how many variants are before and after each gene's start and end coordinate
    start_tally=[]
    end_tally=[]
    for idx,row in excision_df_merged.iterrows():
        cur_gene=row.hgnc_symbol
        cur_vars=common_var_info[cur_gene]
        pre_start_vars=0
        post_end_vars=0
        # if the current gene has some common variants in its excision window...
        if len(cur_vars)>0:
            # check how many are prior to the gene start
            for item in cur_vars:
                if item[0]<row.start_position:
                    pre_start_vars+=1
                elif item[0]>row.end_position:
                    post_end_vars+=1
        start_tally.append(pre_start_vars)
        end_tally.append(post_end_vars)

    excision_df_merged['num_vars_before_start']=start_tally
    excision_df_merged['num_vars_after_end']=end_tally


    # save the information
    excision_df_merged.to_csv(output_dir + '/CommonVars_ALL_summary.txt',sep='\t')
    excision_df_merged.to_csv(output_dir + "/CommonVars_ALL_summary_noIDX.txt", sep='\t', header=None, index=None)
    with open(output_dir + '/CommonVars_ALL_dict.pkl','wb') as file:
        pickle.dump(common_var_info,file)


    # code that was aiming to remove snps preemptively, but we'll do this in later stages of the pipeline
    # new code----------
    # remove any variants that will never encompass a gene exon
    # if (nearby_gene_filter=='True') or (nearby_gene_filter==True) or (nearby_gene_filter==1):
    #     overall_cv_dict={}
    #     for gene in excision_df_merged.hgnc_symbol.unique():
    #         # query the common vars passing the af threshold by gene
    #         pos_list=list(pd.DataFrame(common_var_info[gene],columns=['pos','af'])['pos'])
    #         af_list=list(pd.DataFrame(common_var_info[gene],columns=['pos','af'])['af'])
            
    #         # get all possible snp combinations
    #         snp_pairings=[(min(a, b), max(a, b)) for a, b in itertools.combinations(pos_list, 2)]

    #         # get possible exonic regions across all transcripts of the gene (so, this is not considering vars that occur in every single transcript)
    #         cur_gene_exons=exon_df[exon_df['hgnc_symbol']==gene]
    #         exon_ranges=[]
    #         for idx,row in cur_gene_exons.iterrows():
    #             exon_ranges.append((row.exon_chrom_start,row.exon_chrom_end))
    #         # get the unique ranges
    #         exon_ranges=list(set(exon_ranges))
    #         # sort these ranges by the first value of each tuple
    #         exon_ranges.sort(key=lambda x: min(x[0], x[1]))

    #         # loop through each pair to see if it encompasses at least one exon
    #         pair_encompasses_exon=[]
    #         for pair in snp_pairings:
    #             snps_encompass_exon=False
    #             for e_range in exon_ranges:
    #                 if encompasses(e_range,pair):
    #                     snps_encompass_exon=True
    #             pair_encompasses_exon.append(snps_encompass_exon)

    #         # combine into a pair-boolean dataframe
    #         by_pair_df = pd.DataFrame({
    #             'snp_pair':snp_pairings,
    #             'pair_encompasses_exon':pair_encompasses_exon
    #         })
    #         # get the snps pairs that passed the filter
    #         by_pair_filt = by_pair_df[by_pair_df['pair_encompasses_exon']==True]
    #         passing_snps = list(by_pair_filt['snp_pair'])
    #         #unzip
    #         snp1, snp2 = zip(*passing_snps)
    #         # combine and get unique passing snp values
    #         passing_snps_unique = list(set(list(snp1)+list(snp2)))

    #         # now, loop through the first list of positions and, if a snp is not ever part of a pair that encompasses an exon, remove it from consideration
    #         keep_snp=[]
    #         for pos in pos_list:
    #             if pos in passing_snps_unique:
    #                 keep_snp.append(True)
    #             else:
    #                 keep_snp.append(False)
    #         passing_snps_df = pd.DataFrame({
    #             'pos':pos_list,
    #             'af':af_list,
    #             'keep_snp':keep_snp
    #         })
    #         passing_snps_df = passing_snps_df[passing_snps_df['keep_snp']==True]

    #         # now, store the passing snps in a new dictionary
    #         cur_cv_list=[]
    #         for idx, row in passing_snps_df.iterrows():
    #             cur_cv_list.append([row.pos,row.af])
    #         overall_cv_dict[gene]=cur_cv_list

    #     # see how many variants were lost in this dictionary compared to the original dictionary
    #     num_old_vars=[]
    #     num_refined_vars=[]
    #     var_difs=[]
    #     gene_list=[]
    #     for gene in list(overall_cv_dict.keys()):
    #         gene_list.append(gene)
    #         refined_vars = overall_cv_dict[gene]
    #         old_vars = common_var_info[gene]
    #         var_dif = len(old_vars)-len(refined_vars)

    #         num_old_vars.append(len(old_vars))
    #         num_refined_vars.append(len(refined_vars))
    #         var_difs.append(var_dif)

    #     # combine into df
    #     dif_df = pd.DataFrame({
    #         'hgnc_symbol':gene_list,
    #         'num_old_vars':num_old_vars,
    #         'num_refined_vars':num_refined_vars,
    #         'vars_lost':var_difs
    #     })
    #     dif_df['prop_original_vars_lost']=dif_df['vars_lost']/dif_df['num_old_vars']
    #     # save this information
    #     dif_df.to_csv(output_dir + "/vars_not_encompassing_exons.csv")

    #     # create summary file of same structure
    #     refined_merged=excision_df_merged.merge(dif_df, on='hgnc_symbol', how='left')
    #     # remove old information
    #     refined_merged.drop(labels=['num_common_vars_in_excision_window','num_vars_before_start','num_vars_after_end'],axis=1,inplace=True)
    #     # reorder and rename columns to get the refined data frame to match the excision_df_merged structure
    #     refined_merged = refined_merged[['hgnc_symbol','num_refined_vars','chromosome_name','start_position','end_position']]
    #     refined_merged.rename(columns={'num_refined_vars':'num_common_vars_in_excision_window'},inplace=True)

    #     # now calculate how many vars surround the gene completely
    #     start_tally=[]
    #     end_tally=[]
    #     for idx,row in refined_merged.iterrows():
    #         cur_gene=row.hgnc_symbol
    #         cur_vars=overall_cv_dict[cur_gene]
    #         pre_start_vars=0
    #         post_end_vars=0
    #         # if the current gene has some common variants in its excision window...
    #         if len(cur_vars)>0:
    #             # check how many are prior to the gene start
    #             for item in cur_vars:
    #                 if item[0]<row.start_position:
    #                     pre_start_vars+=1
    #                 elif item[0]>row.end_position:
    #                     post_end_vars+=1
    #         start_tally.append(pre_start_vars)
    #         end_tally.append(post_end_vars)

    #     refined_merged['num_vars_before_start']=start_tally
    #     refined_merged['num_vars_after_end']=end_tally

    #     # save the information
    #     refined_merged.to_csv(output_dir + '/CommonVars_ALL_summary.txt',sep='\t')
    #     refined_merged.to_csv(output_dir + "/CommonVars_ALL_summary_noIDX.txt", sep='\t', header=None, index=None)
    #     with open(output_dir + '/CommonVars_ALL_dict.pkl','wb') as file:
    #         pickle.dump(overall_cv_dict,file)
        

    # # new code above--------







if __name__=='__main__':
    main()

