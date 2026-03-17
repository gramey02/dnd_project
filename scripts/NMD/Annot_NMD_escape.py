### Script to annotate on a per-transcript basis which indels will induce non-sense mediated decay
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle
from itertools import product

# parse any arguments
def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--exon_file', type = str, required = True, help = 'File that contains gene, transcript, and exon information.')
    parser.add_argument('--targetable_common_var_file', type = str, required = True, help = 'File that contains targetable common variant positions by gene.')
    parser.add_argument('--penultimate_rule', type = int, required = True, help = 'Number of bases to consider at end of penultimate exon for escaping NMD.')
    parser.add_argument('--cds_rule', type = int, required = True, help = 'Number of bases to consider downstream of coding start site for escaping NMD.')
    parser.add_argument('--exon_length_rule', type=int, required=True, help='Exon length rule--if a var falls in an exon longer than this, NMD will not be induced.')
    parser.add_argument('--output_dir', type=str,required=True, help='Output directory.')
    args = parser.parse_args()
    return args

# main
def main():
    args = parse_args() # get the arguments passed into the shell script
    # parse the input arguments to the python script here
    exon_file=args.exon_file
    targetable_common_var_file=args.targetable_common_var_file
    penultimate_rule=args.penultimate_rule
    cds_rule=args.cds_rule
    exon_length_rule=args.exon_length_rule
    output_dir=args.output_dir

    # note the NMD escape rules we have to be aware of
    # 1. if the variant falls in the last exon
    # 2. if the variant falls in the last X bases of the penultimate exon (usually 50bp)
    # 3. if the variant falls in the first Y bases of the coding start site (usually 150bp)
    # 4. if the transcript has a single exon
    # 5. if the variant is an exon more than a certain length (usually 407bp)
    # if none of these are satisfied, the transcript should still undergo NMD

    # load data
    exon_df = pd.read_csv(exon_file,index_col=0,dtype={'chromosome_name':str})

    # load and reformat common vars
    with open(targetable_common_var_file, 'rb') as fp:
        targetable_common_var_dict=pickle.load(fp)
    cd_vars_dict={}
    for gene, var_dict in targetable_common_var_dict.items():
        cur_var_pos=[]
        if len(var_dict)>0:
            if type(var_dict[0])!=type('string'):
                for var in var_dict:
                    cur_var_pos.append(var[0])
                cd_vars_dict[gene]=cur_var_pos

    # data preprocessing
    protein_coding_biotypes=['protein_coding','nonsense_mediated_decay', 'non_stop_decay', 'lncRNA', 'miRNA']
    exon_df=exon_df[exon_df.transcript_biotype.isin(protein_coding_biotypes)].reset_index(drop=True) # ge_cds = gene_exons, coding sequences

    # loop over each gene to see if it has variants that can/cannot induce NMD
    gene_dict={}
    for gene in exon_df.hgnc_symbol.unique():
        if gene in list(cd_vars_dict.keys()):
            gene_df = exon_df[exon_df.hgnc_symbol==gene]
            # get the common variant positions for this gene
            cur_cd_vars = cd_vars_dict[gene]
            escapes_NMD=False # this is a variable that OVERALL will represent if a transcript can have NMD induced with at least one of the common exonic coding vars
            # create dictionary to hold transcript information
            transcript_dict={}
            for transcript in gene_df.ensembl_transcript_id.unique():
                # get the current transcript's exon information
                transcript_df = gene_df[gene_df.ensembl_transcript_id==transcript]
                total_exons = max(transcript_df['rank']) # this gets the "#" of the last exon, telling us how many exons there are
                cur_strand = transcript_df.strand.values[0]
                # create variables to store data
                satisfies_1 = []
                satisfies_2 = []
                satisfies_3 = []
                satisfies_4 = []
                satisfies_5 = []
                induces_NMD = []

                # perhaps counterintuitively, let's start with rule #4, since the other rules depend on there being more than one exon
                if total_exons<=1:
                    escapes_NMD=True
                    for var in cur_cd_vars:
                        satisfies_4.append(var)
                else:
                    # check rule #1 - if any variants fall in the last exon
                    # ---------------------------------------------------
                    last_exon_start = (transcript_df[transcript_df['rank']==total_exons])['exon_chrom_start'].values[0] # get the range of this last exon
                    last_exon_end = (transcript_df[transcript_df['rank']==total_exons])['exon_chrom_end'].values[0]
                    # check if common vars fall in this range
                    for var in cur_cd_vars:
                        pos = var
                        if (pos>=last_exon_start and pos<=last_exon_end):
                            escapes_NMD=True
                            satisfies_1.append(pos)
                    
                    # check rule #2 - if any variants fall in the last X bases of the penultimate exon (stl='second to last')
                    # -----------------------------------------------------
                    stl_exon_start = (transcript_df[transcript_df['rank']==total_exons-1])['exon_chrom_start'].values[0] # get the second to last exon coordinates
                    stl_exon_end = (transcript_df[transcript_df['rank']==total_exons-1])['exon_chrom_end'].values[0]
                    if cur_strand==1:
                        lower_bound_to_check = stl_exon_end - penultimate_rule
                        upper_bound_to_check = stl_exon_end
                    else:
                        lower_bound_to_check = stl_exon_start
                        upper_bound_to_check = stl_exon_start + penultimate_rule
                    for var in cur_cd_vars:
                        pos = var
                        if (pos>=lower_bound_to_check and pos<=upper_bound_to_check):
                            escapes_NMD=True
                            satisfies_2.append(pos)
                    
                    # check rule #3 - if any variants fall in the first Y bases after the coding start site
                    # ----------------------------------------------------
                    if cur_strand==1:
                        data = (transcript_df.genomic_coding_start.values)
                        data = data[~np.isnan(data)]
                        genomic_cdss = min(data)
                        lower_bound_to_check = genomic_cdss
                        upper_bound_to_check = genomic_cdss + cds_rule
                    else:
                        data = (transcript_df.genomic_coding_end.values)
                        data = data[~np.isnan(data)]
                        genomic_cdss = max(data)
                        lower_bound_to_check = genomic_cdss - cds_rule
                        upper_bound_to_check = genomic_cdss
                    for var in cur_cd_vars:
                        pos=var
                        if (pos>=lower_bound_to_check and pos<=upper_bound_to_check):
                            escapes_NMD=True
                            satisfies_3.append(pos)
                    
                    # check rule #5 - if the variant falls in an exon larger than a certain amount (default 407bp), NMD is escaped
                    # -----------------------------------------------------
                    # get the exon each var falls in for the current transcript
                    for var in cur_cd_vars:
                        pos=var
                        exon_num=(transcript_df.loc[(transcript_df.exon_chrom_start <= pos) & (pos <= transcript_df.exon_chrom_end), "rank"]).values[0]
                        var_containing_exon_df=transcript_df[transcript_df['rank']==exon_num]
                        var_containing_exon_length=var_containing_exon_df.exon_chrom_end.values[0] - var_containing_exon_df.exon_chrom_start.values[0]
                        if var_containing_exon_length > exon_length_rule:
                            satisfies_5.append(pos)

                    # if the var appears in none of the rule lists, let's say it induces NMD
                    for var in cur_cd_vars:
                        pos=var
                        if ((pos not in satisfies_1) and (pos not in satisfies_2) and (pos not in satisfies_3) and (pos not in satisfies_4) and (pos not in satisfies_5)):
                            induces_NMD.append(pos)

                # calculate some summary statistics to store
                escapes_nmd_set=set()
                escapes_nmd_set.update(satisfies_1, satisfies_2, satisfies_3, satisfies_4, satisfies_5)
                induces_nmd_set=set()
                induces_nmd_set.update(induces_NMD)
                num_vars_considered = len(cur_cd_vars)
                vars_escaping_NMD = len(escapes_nmd_set)
                vars_inducing_NMD = len(induces_nmd_set)

                # store data in transcript dict
                transcript_dict[transcript] = {'num_vars_considered':num_vars_considered,
                                            'num_vars_inducing_NMD':vars_inducing_NMD,
                                            'num_vars_escaping_NMD':vars_escaping_NMD,
                                            'vars_inducing_NMD':induces_nmd_set,
                                            'vars_escaping_NMD':escapes_nmd_set
                }
            # store data in larger gene dict
            gene_dict[gene] = transcript_dict
            # save gene_dict checkpoint
            with open(output_dir + '/NMD/NMD_induction_across_transcripts.pkl','wb') as file:
                pickle.dump(gene_dict,file)

    
    # Now, let's note which vars are have the same effect (NMD/NMD escape) across all transcripts of interest in our set
    # set up some variables
    gene_list=[]
    vars_consistently_inducing_NMD=[]
    vars_consistently_escaping_NMD=[]
    num_consistently_inducing_NMD=[]
    num_consistently_escaping_NMD=[]
    num_considered = []
    for gene,transcript_dict in gene_dict.items():
        # get the total number of transcripts considered for the gene
        total_transcripts=len(transcript_dict.keys())
        # variables storing data
        consistently_induces_NMD=[]
        consistently_escapes_NMD=[]
        # get the vars for the current gene again
        cur_cd_vars = cd_vars_dict[gene]
        if (cur_cd_vars!='No shared exonic regions' and len(cur_cd_vars)!=0):
            for var in cur_cd_vars:
                pos=var
                times_var_induces_NMD=0
                times_var_escapes_NMD=0
                for transcript,summary_dict in transcript_dict.items():
                    if type(summary_dict)==type('string'):
                        #print(summary_dict)
                        continue
                    else:
                        if pos in summary_dict['vars_inducing_NMD']:
                            times_var_induces_NMD+=1
                        if pos in summary_dict['vars_escaping_NMD']:
                            times_var_escapes_NMD+=1
                if times_var_induces_NMD==total_transcripts:
                    consistently_induces_NMD.append(pos)
                if times_var_escapes_NMD==total_transcripts:
                    consistently_escapes_NMD.append(pos)
            num_considered.append(len(cur_cd_vars))
        else:
            num_considered.append(0)
        gene_list.append(gene)
        vars_consistently_inducing_NMD.append(consistently_induces_NMD)
        vars_consistently_escaping_NMD.append(consistently_escapes_NMD)
        num_consistently_inducing_NMD.append(len(consistently_induces_NMD))
        num_consistently_escaping_NMD.append(len(consistently_escapes_NMD))

    # combine into summary data frame
    NMD_summary_df = pd.DataFrame({
        'gene':gene_list,
        'num_vars_considered':num_considered,
        'num_vars_inducing_NMD':num_consistently_inducing_NMD,
        'num_vars_escaping_NMD':num_consistently_escaping_NMD
    })
    # save this data frame
    NMD_summary_df.to_csv(output_dir + '/NMD/NMD_induction_summary.csv',index=None)

    # combine all variant info into a data frame too
    NMD_var_info = pd.DataFrame({
        'gene':gene_list,
        'vars_consistently_inducing_NMD':vars_consistently_inducing_NMD,
        'vars_consistently_escaping_NMD':vars_consistently_escaping_NMD
    })
    NMD_var_info.to_csv(output_dir + '/NMD/NMD_induction_var_info.csv',index=None)

    # filter to only genes that have vars that induce NMD
    genes_to_continue_with = list(NMD_summary_df[NMD_summary_df['num_vars_inducing_NMD']>0]['gene'])
    # combine into dictionary
    pos_dict={}
    NMD_var_info_filt = NMD_var_info[NMD_var_info['gene'].isin(genes_to_continue_with)]
    for idx,row in NMD_var_info_filt.iterrows():
        pos_dict[row.gene]=row.vars_consistently_inducing_NMD

    # create position files for vcf filtering below----------------

    # load vcf files
    vcf_dict={}
    chroms = exon_df.chromosome_name.unique()
    for chrom in chroms:
        af_filename = '/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/TGP_chr' + chrom + '_afs.txt'
        cur_chrom_TGP_afs = pd.read_csv(af_filename, sep=' ', names = ['chrom', 'pos', 'ref', 'alt', 'ac', 'an', 'af', 'afr_af', 'amr_af', 'eas_af', 'eur_af', 'sas_af'])
        vcf_dict[chrom] = cur_chrom_TGP_afs[['chrom','pos','ref','alt']]

    # for each gene, format its common vars where each row contains chrom & pos. We'll filter larger vcfs by this information later using bcftools
    savedir=output_dir + "/excavate/CommonVar_locs/"
    genes_w_commonVars=[]
    lowest_var_pos=[]
    highest_var_pos=[]
    for gene,cv_list in pos_dict.items():
        if (len(cv_list)>0):
            if (type(cv_list[0])!=type('string')):
                # convert the lists into a data frame
                cv_df = pd.DataFrame(cv_list, columns=['pos'])
                cur_chrom = ((exon_df[exon_df['hgnc_symbol']==gene])['chromosome_name'].values[0])
                cv_df['chrom'] = cur_chrom
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

#----------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()