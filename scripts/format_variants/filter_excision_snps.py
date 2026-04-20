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
    parser.add_argument('--cv_dict_filepath', type = str, required = True, help = 'Filepath for the common var dict.')
    parser.add_argument('--exon_file', type = str, required = True, help = 'Exon information file name.')
    parser.add_argument('--output_dir', type=str, required=True, help="output directory for the pipeline run.")
    parser.add_argument('--gene', type=str, required=True, help="Current gene's name.")
    parser.add_argument('--valid_snp_dir', type=str, required=True, help="Location to save valid snp pairs to.")
    args = parser.parse_args()
    return args

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
        if cur_snp_info==False:
            return False
        else:
            transcripts_where_exon_is_excised+=1
    # # if the snp range always encompasses an exon, even across all transcripts, return true; else, return False
    # if ((transcripts_where_exon_is_excised<total_transcripts) and (transcripts_where_exon_is_excised>0)):
    #     print(snp_range)
    #     print(transcript)
    return (transcripts_where_exon_is_excised==total_transcripts)


def main():
    # parse input arguments
    args = parse_args()
    cv_dict_filepath=args.cv_dict_filepath
    exon_file=args.exon_file
    output_dir=args.output_dir
    gene=args.gene
    valid_snp_dir=args.valid_snp_dir

    os.makedirs(output_dir, exist_ok=True)

    # load exon data frame
    exon_df = pd.read_csv(exon_file, dtype={'chromosome_name':'str'},index_col=0)

    # load common var dict
    with open(cv_dict_filepath, 'rb') as fp:
        cv_dict = pickle.load(fp)

    gene_exons=exon_df[exon_df['hgnc_symbol']==gene]
    gene_exons = gene_exons[gene_exons['transcript_biotype'].isin(['protein_coding','nonsense_mediated_decay', 'non_stop_decay', 'lncRNA', 'miRNA'])]

    snp_positions=list(pd.DataFrame(cv_dict[gene], columns=['pos', 'af'])['pos'])
    possible_snp_pairs = set(combinations(snp_positions, 2)) # list(combinations(snp_positions, 2))
    refined_snp_list={x for x in possible_snp_pairs if pair_encompasses_exon_ubiquitously(x,gene_exons)}
    # save the valid snp pairs
    with open(os.path.join(valid_snp_dir, f'{gene}_valid_snp_pairs.pkl'), 'wb') as fp:
        pickle.dump(refined_snp_list, fp)
    # un-pair the snps and save as a final set
    final_snp_set=set()
    for pair in refined_snp_list:
        if pair[0] not in final_snp_set:
            final_snp_set.add(pair[0])
        if pair[1] not in final_snp_set:
            final_snp_set.add(pair[1])
    #final_snp_list=list(final_snp_set)
    # save the new list of snps
    output_filepath = os.path.join(output_dir, f'{gene}_refined_snp_list.pkl')
    with open(output_filepath, 'wb') as fp:
        pickle.dump(final_snp_set,fp)



if __name__=='__main__':
    main()
