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
    # set some variables
    savedir=output_dir + "/excavate/CommonVar_locs/"
    genes_w_commonVars=[]
    lowest_var_pos=[]
    highest_var_pos=[]

    # if the current strategy is excision, run the loop below, taking into account the filtered snps that form valid excision pairs
    if (edit_strat is not None) and (edit_strat=='excision'):
        # load refined common variants
        refined_vars_gene_fps=os.listdir(cv_dict_filepath)
        for file in refined_vars_gene_fps:
            cur_gene=file.split('_')[0]
            with open(os.path.join(cv_dict_filepath, file), 'rb') as fp:
                cur_snps=pickle.load(fp)
            if len(cur_snps)>0:
                cur_chrom = ((exon_df[exon_df['hgnc_symbol']==cur_gene])['chromosome_name'].values[0])
                cv_df = pd.DataFrame(cur_snps, columns=['pos'])
                genes_w_commonVars.append(cur_gene)
                cv_df['chrom']=cur_chrom
                cv_df = cv_df[['chrom','pos']]
                cv_df.to_csv(savedir + cur_gene + '_CommonVar_locs.txt',sep='\t',index=False, header=False)
                # save the upper and lower bounds of where the variants occur, plus some padding to provide ample search space for PAMs
                lowest_var_pos.append(min(cv_df.pos)-100)
                highest_var_pos.append(max(cv_df.pos)+100)

    else:
        # load common var dict
        with open(cv_dict_filepath, 'rb') as fp:
            cv_dict = pickle.load(fp)
        
        # load vcf files
        vcf_dict={}
        chroms = exon_df.chromosome_name.unique()
        for chrom in chroms:
            af_filename = os.path.join(af_file_dir, 'TGP_chr' + chrom + '_afs.txt')
            cur_chrom_TGP_afs = pd.read_csv(af_filename, sep=' ', names = ['chrom', 'pos', 'ref', 'alt', 'ac', 'an', 'af', 'afr_af', 'amr_af', 'eas_af', 'eur_af', 'sas_af'])
            vcf_dict[chrom] = cur_chrom_TGP_afs[['chrom','pos','ref','alt']]
        
        # for each gene, format its common vars where each row contains chrom & pos. We'll filter larger vcfs by this information later using bcftools
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
        'hgnc_symbol': pd.Series(genes_w_commonVars, dtype='object'),
        'lower_coord': lowest_var_pos,
        'higher_coord': highest_var_pos
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
