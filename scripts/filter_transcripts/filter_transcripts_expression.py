import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle

#parse the arguments passed into the script
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--transcript_tpm_file', type = str, required = True, help = 'Filepath to transcript tpms from GTEx.')
    parser.add_argument('--sample_attributes_file', type = str, required = True, help = 'Filepath to attributes of samples in transcript tpm file (from GTEx)')
    parser.add_argument('--gene_median_tpms_file', type = str, required = True, help = 'Filepath to gene tpms "median-ed" across samples from GTEx.')
    parser.add_argument('--exon_file', type = str, required = True, help = 'Exon information file name.')
    parser.add_argument('--tpm_thresh', type = float, required = True, help = 'Transcripts per million threshold to consider something "expressed" in a certain tissue.')
    parser.add_argument('--prop_thresh', type = float, required = True, help = 'Threshold for proportion of expression a transcript must account for to be considered for furuter analyses.')
    parser.add_argument('--keep_all_transcripts', type = str, required = True, help = 'Keep all transcripts for future analyses if no transcripts pass the proportion threshold.')
    parser.add_argument('--output_file', type = str, required = True, help = 'Output file holding filtered data.')
    args = parser.parse_args()
    return args

def main():
    # parse args from shell script
    args = parse_args() # get the arguments passed into the shell script
    gene_tpm_filename = args.gene_median_tpms_file # '/wynton/group/capra/data/GTEx/v10/2025-05-06/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct'
    transcript_tpm_filename = args.transcript_tpm_file # '/wynton/protected/home/capra/gramey02/ConklinCollab/data/GTEx/transcripts_tpm_small.txt'
    exon_filename = args.exon_file
    sample_attributes_filename = args.sample_attributes_file # '/wynton/group/capra/data/GTEx/v10/2025-05-06/metadata/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt'
    output_file=args.output_file

    #load file that contains column names for ttpm
    ttpm_colnames_file="/wynton/protected/home/capra/gramey02/ConklinCollab/data/GTEx/transcripts_tpm_small.txt"
    ttpm_colnames_df=pd.read_table(ttpm_colnames_file)
    # load transcript tpms
    ttpms = pd.read_table(transcript_tpm_filename,names=list(ttpm_colnames_df.columns))
    # load sample attributes
    dtype_dict_small={'SMGTC':'str'}
    sa = pd.read_table(sample_attributes_filename,dtype=dtype_dict_small)

    # transcript tpm file is large, so split into batches respecting gene boundaries
    # list of lists with gene names to look into each time

    # create a dictionary mapping sample id to tissue
    sample_tissue_dict = dict(zip(list(sa.SAMPID), list(sa.SMTSD)))
    # get transcript-->gene mappings
    g_t_map = ttpms[['transcript_id', 'gene_id']]

    # reformat tpm df
    ttpms_t = ttpms.transpose()
    ttpms_t.columns = ttpms_t.iloc[0]
    ttpms_t = ttpms_t[2:]
    ttpms_t.reset_index(inplace=True)

    # map sample ids to tissue names
    ttpms_t = ttpms_t.assign(tissue_names=ttpms_t['index'].map(sample_tissue_dict))
    # reformat before grouping
    ttpms_vals = ttpms_t.drop(labels=['index'], axis=1)
    # group by median value of tpms for each tissue
    tpm_grouped = ttpms_vals.groupby('tissue_names').median()
    # transpose df again
    grouped_t = tpm_grouped.transpose() # now the columns should be tissues, and the rows should be transcripts
    grouped_t.reset_index(inplace=True)
    # add the gene column back in
    gene_transcript_dict = dict(zip(list(g_t_map.transcript_id),list(g_t_map.gene_id)))
    ttpm_byTissue = grouped_t.assign(gene_id=grouped_t['transcript_id'].map(gene_transcript_dict))
    # remove cell lines
    ttpm_byTissue.drop(labels=['Cells - EBV-transformed lymphocytes','Cells - Cultured fibroblasts'],axis=1,inplace=True)
    # ----------------------------------Data Frame is formatted now!----------------------------------

    # load gene expression file
    gene_tpm = pd.read_csv(gene_tpm_filename, sep='\t', skiprows=[0,1])
    # now get the proportion of expression each transcript contributes to each tissue
    tissue_list = [x for x in list(ttpm_byTissue.columns) if (('transcript' not in x) and ('gene' not in x))]
    # run a loop to note the proportion of expression that each transcript contributes to within a tissue
    prop_df = None
    for gene in ttpm_byTissue.gene_id.unique():
        ttpm_filt = ttpm_byTissue[ttpm_byTissue.gene_id==gene]
        tissue_expr_totals = (pd.DataFrame((ttpm_filt.drop(labels=['transcript_id','gene_id'],axis=1)).sum(axis=0), columns=['total_tpm'])).reset_index()
        cur_gene_transcript_props=[]
        cur_gene_transcript_tissueIDs=[]
        for tissue in tissue_list:
            cur_gene_transcript_tissueIDs.append(tissue)
            cur_total = tissue_expr_totals[tissue_expr_totals.tissue_names==tissue]['total_tpm'].values[0]
            if cur_total!=0:
                cur_tissue_props = ttpm_filt[tissue]/cur_total
            else:
                cur_tissue_props = np.zeros((len(ttpm_filt[tissue]),1))
            intermediate_df = pd.DataFrame({
                'tissue':tissue,
                'gene': gene,
                'transcript_id': ttpm_filt.transcript_id,
                'prop_of_expr_within_gene':list(cur_tissue_props)
            })
            prop_df = pd.concat([prop_df,intermediate_df])
    # reset the index of the df
    prop_df.reset_index(inplace=True, drop=True)
    # we'll filter on a specific proportion threshold below
    # -------------------------------------------------------------------------------------------------

    # Now let's determine which of these transcripts is not 'important enough' to include in future analyses, 
    # based on expression

    # load thresholds
    tpm_thresh = args.tpm_thresh
    prop_thresh = args.prop_thresh
    # load ensembl gene information to match up/fill in with GTEx info
    dtype_dict_small={'chromosome_name':'str','start_position':'int'}
    exon_df = pd.read_csv(exon_filename, index_col=0,dtype=dtype_dict_small)
    # parse boolean
    keep_all_transcripts=args.keep_all_transcripts

    # note where each gene is expressed
    transcripts_to_keep=[]
    prop_transcripts_lost_per_gene = []
    tissue_mappings=pd.read_csv('/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/transcript_tpm/tissue_mappings.csv',index_col=0)
    tissue_mapping_dict=dict(zip(tissue_mappings['gene_underscored_tissue'],tissue_mappings['transcript_space_tissue']))
    for gene in prop_df.gene.unique():
        noV_gene = gene.split('.')[0] # get the gene name without the decimal version attached
        # filter the transcript tpm df to just this tissue and gene
        prop_filt = prop_df[prop_df.gene==gene]
        # reformat the gene expression df
        gene_expr_df = (gene_tpm[gene_tpm.Name==gene]).transpose().reset_index()
        gene_expr_df.columns=gene_expr_df.iloc[0]
        gene_expr_df = gene_expr_df[2:]
        gene_expr_df = gene_expr_df[~gene_expr_df.Name.isin(['Cells_Cultured_fibroblasts', 'Cells_EBV-transformed_lymphocytes'])]
        # filter to just tissues the gene is expressed in
        gene_expr_tissues = gene_expr_df[gene_expr_df[gene]>=tpm_thresh]
        expr_tissue_list = list(gene_expr_tissues.Name)
        mapped_tissues=[]
        for val in expr_tissue_list:
            mapped_tissues.append(tissue_mapping_dict[val])
        # filter the prop_df to just these tissues the gene is expressed in
        prop_expr_tissues = prop_filt[prop_filt.tissue.isin(mapped_tissues)]
        # now get each transcript's median expression level across tissues the gene is expressed in
        prop_grouped = (pd.DataFrame(prop_expr_tissues.groupby('transcript_id')['prop_of_expr_within_gene'].mean())).reset_index() # median())).reset_index()
        total_transcripts = len(prop_grouped.transcript_id)
        prop_grouped_below_thresh = prop_grouped[prop_grouped['prop_of_expr_within_gene']<prop_thresh]
        if total_transcripts>0:
            prop_transcripts_not_passing_thresh = len(prop_grouped_below_thresh.transcript_id)/total_transcripts
            # save this ratio
            prop_transcripts_lost_per_gene.append(prop_transcripts_not_passing_thresh)
        # get the transcripts that do pass the threshold
        prop_grouped_above_thresh = prop_grouped[prop_grouped['prop_of_expr_within_gene']>=prop_thresh]
        if len(prop_grouped_above_thresh.transcript_id.unique())==0:
            exon_filt = exon_df[exon_df.ensembl_gene_id==noV_gene]
            if keep_all_transcripts=="True":
                # keep all transcripts for the gene listed in ensembl
                for transcript in exon_filt.ensembl_transcript_id.unique():
                    transcripts_to_keep.append(transcript)
            else:
                # keep only the canonical transcript if "keep_all" flag is false
                canonical = exon_filt[exon_filt.transcript_is_canonical==1]
                for transcript in canonical.ensembl_transcript_id.unique():
                    transcripts_to_keep.append(transcript)
        else:
            for transcript in prop_grouped_above_thresh.transcript_id.unique():
                transcripts_to_keep.append(transcript)
    #--------------------------------------------------------------------------------------------------

    # Now we have a list of transcripts that we want to move forward with in future analyses.
    # Filter the exon dataframe by these transcripts and pass into future analyses.

    # first need to separate out the '.version' from some of the transcript names
    noV_transcripts = []
    for transcript in transcripts_to_keep:
        if '.' in transcript:
            noV_transcripts.append(transcript.split('.')[0])
        else:
            noV_transcripts.append(transcript)
    # now filter the exon df to these transcripts only
    final_exon_df = exon_df[exon_df.ensembl_transcript_id.isin(noV_transcripts)]

    # append all transcripts to the data frame if the gene is completely missing
    missing_genes = [x for x in list(exon_df.ensembl_gene_id.unique()) if x not in list(final_exon_df.ensembl_gene_id.unique())]
    for gene in missing_genes:
        noV_transcripts = noV_transcripts+ list(exon_df[exon_df.ensembl_gene_id==gene]['ensembl_transcript_id'].unique())
    # refilter the exon df after adding these transcripts
    final_exon_df = exon_df[exon_df.ensembl_transcript_id.isin(noV_transcripts)]

    # save file
    print(output_file)
    final_exon_df.to_csv(output_file)

#-----------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()

