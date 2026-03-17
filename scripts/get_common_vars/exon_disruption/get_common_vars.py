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
    parser.add_argument('--af_limit', type = float, required = True, help = 'Lower bound of allele frequency for common variants.')
    parser.add_argument('--chrom', type = str, required = True, help = 'Chromosome number.')
    parser.add_argument('--af_file', type = str, required = True, help = 'Allele frequency file name (file contains allele frequencies across populations by position).')
    parser.add_argument('--exon_file', type = str, required = True, help = 'Exon information file name.')
    parser.add_argument('--output_dir', type=str, required=True, help="output directory for the pipeline run.")
    parser.add_argument('--total_num_chroms', type=int, required=True, help="Total number of unique chromosomes for gene set of interest.")
    args = parser.parse_args()
    return args

def main():
    args = parse_args() # get the arguments passed into the shell script
    # parse the input arguments to the python script here
    af_limit=args.af_limit
    chrom=args.chrom
    af_file=args.af_file
    exon_file=args.exon_file
    shared_exon_savepath=args.output_dir + '/ubiq_regions/'
    common_vars_savepath=args.output_dir + '/ubiq_region_CommonVars/'
    total_num_chroms=args.total_num_chroms

    # load allele frequencies for snps around the gene
    snp_afs=pd.read_table(af_file, sep=' ', names = ['chr', 'pos', 'ref', 'alt', 'ac', 'an', 'af', 'afr_af', 'amr_af', 'eas_af', 'eur_af', 'sas_af'])
    # filter to those that satisfy the allele frequency threshold
    af_pass=snp_afs[(snp_afs.af>=af_limit) & (snp_afs.af<=1-af_limit)]
    # get a list of common variant positions in and around the gene
    pos=list(af_pass.pos)
    # df of exonic information across genes and transcripts
    dtype_dict_small={'chromosome_name':'str'} # set type of one of the columns, since it's ambiguous atm
    exon_df=pd.read_csv(exon_file,index_col=0,dtype=dtype_dict_small)
    chrom_genes=exon_df[exon_df.chromosome_name==chrom]


    #----new code taking exonic region overlap
    # create dictionary that will hold exonic information on each gene
    exonic_range_byGene={}
    for gene in chrom_genes.hgnc_symbol.unique():
        gene_exons=exon_df[exon_df['hgnc_symbol']==gene]
        protein_coding_biotypes=['protein_coding','nonsense_mediated_decay', 'non_stop_decay', 'lncRNA', 'miRNA']
        ge_cds=gene_exons[gene_exons.transcript_biotype.isin(protein_coding_biotypes)].reset_index(drop=True) # ge_cds = gene_exons, coding sequences
        # run a loop getting the exonic regions of each transcript (later we'll find the overlap among all of them)
        singleGene_transcript_exon_ranges={}
        for transcript in ge_cds['ensembl_transcript_id'].unique():
            transcript_df=ge_cds[ge_cds['ensembl_transcript_id']==transcript] # we can also have a flag here to ignore transcripts that aren't commonly expressed enough
            # get the total number of transcripts for the gene
            num_transcripts=len(transcript_df.ensembl_transcript_id.unique())
            # create a dict to hold the exon ranges for each exon
            cur_transcript_exon_ranges=[]
            for idx,row in transcript_df.iterrows():
                cur_exon=row.ensembl_exon_id
                exon_start = row.exon_chrom_start
                exon_end = row.exon_chrom_end
                # store in dict
                cur_transcript_exon_ranges.append([cur_exon,exon_start,exon_end])
            # store the current exonic ranges among the ranges for the other transcripts (turn it into a dataframe first)
            singleGene_transcript_exon_ranges[transcript]=pd.DataFrame(cur_transcript_exon_ranges, columns=['ensembl_exon_id','exon_start','exon_end'])
        # add to larger dictionary for each gene
        exonic_range_byGene[gene] = singleGene_transcript_exon_ranges



    # now calculate the shared exonic regions across all transcripts per gene
    universal_exonic_regions_byGene={}
    for gene,transcript_dict in exonic_range_byGene.items():
        # get number of transcripts for the current gene
        transcript_list = list(transcript_dict.keys())
        num_transcripts = len(transcript_list)
        # initialize the first range to compare to as the first transcript's exon regions
        if len(transcript_list)>0:
            transcript1 = transcript_list[0]
            t1_df = transcript_dict[transcript_list[0]]
            t1_df['start_end'] = t1_df[['exon_start','exon_end']].apply(tuple,axis=1)
            t1_exon_loc_list = list(t1_df['start_end'])
            # have an if statement here checking if there is only one transcript for the gene (otherwise don't need to run all of the stuff below)
            for transcript2, t2_df in transcript_dict.items():
                if len(t1_exon_loc_list)!=0:
                    if transcript1!=transcript2:
                        # if the transcripts we're comparing aren't the same, compare their exon ranges, one exon of t1 at a time
                        t2_copy = t2_df.copy()
                        t2_num_exons = len(t2_copy.ensembl_exon_id)
                        
                        # turn start and end into tuples
                        t2_copy['start_end'] = t2_copy[['exon_start','exon_end']].apply(tuple,axis=1)
                        t2_exon_loc_list = list(t2_copy['start_end'])
                        
                        # create a product to get every combination of exon ranges that can be compared
                        range_combos = list(product(t2_exon_loc_list, t1_exon_loc_list))
                        combos_df = pd.DataFrame(range_combos,columns=['t2_vals','t1_vals'])
                        combos_df[['t2a', 't2b']] = combos_df.t2_vals.tolist()
                        combos_df[['t1a','t1b']] = combos_df.t1_vals.tolist()
                        combos_df = combos_df[['t2a','t2b','t1a','t1b']]
                        
                        # calculate overlap between each combo
                        combos_df['overlap_start']=combos_df[['t1a','t2a']].max(axis=1)
                        combos_df['overlap_end']=combos_df[['t1b','t2b']].min(axis=1)
                        combos_df['overlap'] = combos_df.overlap_end - combos_df.overlap_start
                        combos_df['overlap_bool'] = combos_df.overlap>0 # we'll consider 0 overlapping because this means a single bp is shared
                        
                        # now let's filter to just the regions that overlap
                        o_regions = combos_df[combos_df.overlap_bool==True]
                        # convert overlapping regions back into tuples
                        o_regions = o_regions.assign(overlapping_regions = o_regions[['overlap_start','overlap_end']].apply(tuple,axis=1))
        
                        # now reset t1_exon_loc_list to compare to the rest of the transcripts
                        t1_exon_loc_list = list(o_regions['overlapping_regions'])
                else:
                    break
        else:
            t1_exon_loc_list=[]
                
        # save the 'universal exonic regions' that we found
        universal_exonic_regions_byGene[gene] = t1_exon_loc_list
        with open(shared_exon_savepath + 'ubiq_exons_chr' + chrom + '.pkl','wb') as file:
            pickle.dump(universal_exonic_regions_byGene, file)
    
    #-----------------------------------------------------------------------------------------
            
    # now that we have the shared exonic regions, we can see what common variants fall in those regions
    common_var_info={}
    num_common_vars_in_exons=[]
    genes=[]
    for gene,exonic_regions in universal_exonic_regions_byGene.items():
        genes.append(gene)
        cur_gene_common_var_info=[]
        if len(exonic_regions)!=0:
            # turn list into a data frame
            for region in exonic_regions:
                region_start=region[0]
                region_end=region[1]
                in_exonic_regions=[x for x in pos if (x>=region_start and x<=region_end)]
                # get the allele frequency information for these positions too
                pos_afs = af_pass[af_pass['pos'].isin(in_exonic_regions)][['pos','af']]
                for idx,row in pos_afs.iterrows():
                    cur_gene_common_var_info.append([row.pos,row.af])
            num_common_vars_in_exons.append(len(cur_gene_common_var_info))
        else:
            num_common_vars_in_exons.append(0)
            cur_gene_common_var_info = ['No shared exonic regions']
        common_var_info[gene] = cur_gene_common_var_info
    # combine summary information into data frame
    exonic_cv_df = pd.DataFrame({
        'gene':genes,
        'num_common_vars_in_shared_exons':num_common_vars_in_exons
    })
    # save the information here
    exonic_cv_df.to_csv(common_vars_savepath + 'CommonVars_chr' + chrom + '_summary.txt',sep='\t')
    with open(common_vars_savepath + '/CommonVars_chr' + chrom + '_dict.pkl','wb') as file:
        pickle.dump(common_var_info,file)


    # merge separate chromosome files
    if len(os.listdir(common_vars_savepath))==(total_num_chroms*2):
        directory=common_vars_savepath
        summary_files=os.listdir(directory)
        summary_files=[x for x in summary_files if 'summary' in x]
        dict_files=os.listdir(directory)
        dict_files=[x for x in dict_files if 'dict' in x]
        summary_df=None
        for file in summary_files:
            int_string=file.split("CommonVars_chr")[1]
            cur_chrom=int_string.split("_summary.txt")[0]
            cur_df=pd.read_table(directory+'/'+file, index_col=0)
            cur_df['chrom']=cur_chrom
            summary_df = pd.concat([summary_df,cur_df])
        summary_df.to_csv(directory + '/CommonVars_ALL_summary.txt', sep='\t')
        summary_df.to_csv(directory + '/CommonVars_ALL_summary_noIDX.txt', sep='\t', header=False, index=False)
        full_dict={}
        for file in dict_files:
            with open(directory + '/' + file, 'rb') as fp:
                full_dict.update(pickle.load(fp))
        with open(directory+'/CommonVars_ALL_dict.pkl', 'wb') as fp:
            pickle.dump(full_dict,fp)


#----------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
    
