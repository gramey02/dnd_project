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
    parser.add_argument('--guides_filepath', type = str, required = True, help = 'Filepath for the common var dict.')
    parser.add_argument('--exon_file', type = str, required = True, help = 'Exon information file name.')
    parser.add_argument('--output_dir', type=str, required=True, help="output directory for the pipeline run.")
    args = parser.parse_args()
    return args

def main():
    # parse input args
    args = parse_args()
    guides_filepath=args.guides_filepath
    exon_file=args.exon_file
    output_dir=args.output_dir

    # load exon file
    exon_df=pd.read_csv(exon_file, index_col=0)

    # load directory names in the excavate output location
    genes_w_guides=[]
    chroms=[]
    gene_dirs=os.listdir(guides_filepath)
    for item in gene_dirs:
        # check if the directory is empty
        if len(os.listdir(guides_filepath + "/" + item)) > 0:
            # check if the files in the directory are empty
            cur_guides = pd.read_csv(guides_filepath + "/" + item + "/" + "all_guides.csv") # load file
            if len(cur_guides['SNP position']) > 0:
                cur_gene=item.split("_")[0]
                cur_chrom=exon_df[exon_df['hgnc_symbol']==cur_gene]['chromosome_name'].values[0]
                genes_w_guides.append(cur_gene)
                chroms.append(cur_chrom)
                # create a text file with the positions of these guides
                cur_guide_info = pd.DataFrame({
                    'chrom':cur_chrom,
                    'pos':list(cur_guides['SNP position'].unique())
                })
                cur_guide_info.to_csv(output_dir + "/excavate/Guide_locs/" + cur_gene + "_Guide_locs.txt", sep="\t", index=False, header=False)

    # save genes that have valid guides
    valid_guides=pd.DataFrame({
        'gene':genes_w_guides,
        'chrom':chroms
    })
    valid_guides.to_csv(output_dir + "/excavate/het_individuals/metadata/" + "genes_w_valid_guides.txt", sep='\t',header=False, index=False)


if __name__=='__main__':
    main()