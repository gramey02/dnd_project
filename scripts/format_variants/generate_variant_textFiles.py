import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle
from itertools import product

#parse the arguments passed into the script
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--cv_dict_filepath', type = str, required = True, help = 'Filepath for the common var dict.')
    parser.add_argument('--exon_file', type = str, required = True, help = 'Exon information file name.')
    parser.add_argument('--output_dir', type=str, required=True, help="output directory for the pipeline run.")
    parser.add_argument('--af_file_dir', type=str, required=True, help="Directory containing TGP_chr*_afs.txt files.")
    parser.add_argument('--edit_strat', type=str, required=False, help="Optional argument that allows for additional filtering to the variant/gene set based on the current edit strategy.")
    args = parser.parse_args()
    return args

def main():
    # parse input arguments
    args = parse_args()
    cv_dict_filepath=args.cv_dict_filepath
    exon_file=args.exon_file
    output_dir=args.output_dir
    af_file_dir=args.af_file_dir
    edit_strat=args.edit_strat

    # load exon data frame
    exon_df = pd.read_csv(exon_file, dtype={'chromosome_name':'str'},index_col=0)

    # load common var dict
    with open(cv_dict_filepath, 'rb') as fp:
        cv_dict = pickle.load(fp)
    # # if the editing strategy requires extra filtering of the set, do it here
    # if edit_strat=="excision":
    #     exon_filt=exon_df[['hgnc_symbol','chromosome_name','start_position','end_position']].drop_duplicates()
    #     genes_to_remove=[]
    #     for idx,row in exon_filt.iterrows():
    #         if row.hgnc_symbol in cv_dict.keys():
    #             cv_list=cv_dict[row.hgnc_symbol]
    #             num_vars_before_gene_start=0
    #             num_vars_after_gene_end=0
    #             for item in cv_list:
    #                 if item[0]<row.start_position:
    #                     num_vars_before_gene_start+=1
    #                 elif item[0]>row.end_position:
    #                     num_vars_after_gene_end+=1
    #             if (num_vars_before_gene_start<1) or (num_vars_after_gene_end<1):
    #                 genes_to_remove.append(row.hgnc_symbol)
    #     # remove genes that don't have enough variants for excision from the dictionary, so they're not passed through to excavate
    #     for k in genes_to_remove:
    #         cv_dict.pop(k, None)

    
    # load vcf files
    vcf_dict={}
    chroms = exon_df.chromosome_name.unique()
    for chrom in chroms:
        af_filename = os.path.join(af_file_dir, 'TGP_chr' + chrom + '_afs.txt')
        cur_chrom_TGP_afs = pd.read_csv(af_filename, sep=' ', names = ['chrom', 'pos', 'ref', 'alt', 'ac', 'an', 'af', 'afr_af', 'amr_af', 'eas_af', 'eur_af', 'sas_af'])
        vcf_dict[chrom] = cur_chrom_TGP_afs[['chrom','pos','ref','alt']]
    
    # for each gene, format its common vars where each row contains chrom & pos. We'll filter larger vcfs by this information later using bcftools
    savedir=output_dir + "/excavate/CommonVar_locs/"
    genes_w_commonVars=[]
    lowest_var_pos=[]
    highest_var_pos=[]
    for gene,cv_list in cv_dict.items():
        if (len(cv_list)>0):
            if (type(cv_list[0])!=type('string')):
                # convert the lists into a data frame
                cv_df = pd.DataFrame(cv_list, columns=['pos','af'])
                cur_chrom = ((exon_df[exon_df['hgnc_symbol']==gene])['chromosome_name'].values[0])
                cv_df['chrom'] = cur_chrom
                cv_df = cv_df.drop(labels=['af'],axis=1)
                cv_df=cv_df[['chrom','pos']]
                # this df of vars should already be filtered to the previously set af_threshold in the get_common_exonic_vars script
                genes_w_commonVars.append(gene) # save the gene for later filtering
                # save this df now for filtering vcfs
                cv_df.to_csv(savedir + gene + '_CommonVar_locs.txt',sep='\t',index=False, header=False)
                # save the upper and lower bounds of where the variants occur, plus some padding to provide ample search space for PAMs
                lowest_var_pos.append(min(cv_df.pos)-100)
                highest_var_pos.append(max(cv_df.pos)+100)

    # also create a file that contains start and end coordinates of each gene, as this is also a necessary input for excavate
    savedir2=output_dir + "/excavate/input_metadata/"
    exon_filt = exon_df[['hgnc_symbol','chromosome_name']].drop_duplicates().reset_index(drop=True)
    # create a data frame of coordinates
    coord_df = pd.DataFrame({
        'hgnc_symbol':genes_w_commonVars,
        'lower_coord':lowest_var_pos,
        'higher_coord':highest_var_pos
    })
    # filter to only those genes that have a valid common variant that we can run EXCAVATE on
    genes_w_commonVars_df = exon_filt[exon_filt['hgnc_symbol'].isin(genes_w_commonVars)]
    # merge with coordinate data
    merged = genes_w_commonVars_df.merge(coord_df, on='hgnc_symbol', how='left')
    merged['coords']='chr'+ merged.chromosome_name.astype(str) + ":" + merged.lower_coord.astype(int).astype(str) + "-" + merged.higher_coord.astype(int).astype(str)
    # add some additional metadata
    merged['filtering_txt_filepath'] = savedir + merged['hgnc_symbol'] + '_CommonVar_locs.txt'
    merged['excavate_input_vcf_filepath'] = output_dir+"/excavate/input_vcfs/" + merged['hgnc_symbol'] + "_input_vcf.gz"
    # save
    merged.to_csv(savedir2 + 'excavate_run_metadata.txt', sep='\t', header=False, index=False)

#----------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
