#!/usr/bin/env python3
"""Generated from Creating_Master_DnD_DataFrame.ipynb."""

# %% [markdown] cell 1
# # Master Data Frame Creation
# In this notebook, we'll combine many data sources to get information on D&D genes all in one place

# %% cell 2
import argparse
import os
import sys
from pathlib import Path

import networkx as nx
import numpy as np
import obonet
import pandas as pd
from datetime import datetime
from itertools import permutations, combinations, product
import itertools
DATE = datetime.now().strftime("%Y-%m-%d")

import matplotlib.pyplot as plt
#import seaborn as sns
from matplotlib.ticker import FuncFormatter
import pickle
import math

# %% [markdown] cell 3
# ## Set up results filepath

# %% cell 4
SCRIPT_DIR = Path(__file__).resolve().parent
DND_PROJECT_ROOT = SCRIPT_DIR.parent.parent
WORKSPACE_ROOT = DND_PROJECT_ROOT.parent
CONKLIN_ROOT = WORKSPACE_ROOT / "ConklinCollab"
CLINGEN_DIR = DND_PROJECT_ROOT / "data" / "clingen"
DND_HGNC_DIR = DND_PROJECT_ROOT / "data" / "dnd_hgnc"

parser = argparse.ArgumentParser(
    description="Build the master D&D summary dataframe from pipeline outputs."
)
parser.add_argument(
    "--results-dir",
    default=str(CONKLIN_ROOT / "data/dHS_and_related_GeneSets/pipeline_results"),
    help="Base directory containing run outputs.",
)
parser.add_argument(
    "--run-name",
    default="RUN_12_19_25",
    help="Pipeline run name to summarize.",
)
parser.add_argument(
    "--output-csv",
    default=None,
    help="Optional path for the final output CSV. Defaults to summary_files/dnd_593_<DATE>.csv under the selected run.",
)
args = parser.parse_args()

results_dir = args.results_dir
run_name = args.run_name
checkpoint_csv = os.path.join(results_dir, run_name, "summary_files", "dnd_593.csv")
output_csv = args.output_csv or os.path.join(
    results_dir,
    run_name,
    "summary_files",
    f"dnd_593_{DATE}.csv",
)

# %% [markdown] cell 5
# ## Load D&D gene set

# %% cell 6
exon_file=os.path.join(results_dir,run_name,"filtered_transcripts/filtered_exon_info.csv")
exon_df=pd.read_csv(exon_file,index_col=0,dtype={"chromosome_name":"str"})

# %% [markdown] cell 7
# ## Load disease and dosage information

# %% cell 8
gd_fp = CLINGEN_DIR / "Clingen-Gene-Disease-Summary-2025-11-13.csv" # gene-disease table from clingen
ds_fp = CLINGEN_DIR / "Clingen-Dosage-Sensitivity-2025-11-13.csv" # dosage sensitivity table from clingen

# load gd table
gd = pd.read_csv(gd_fp)
cols = gd.iloc[3].values
gd = gd[5:]
gd.columns = cols

# load dosage sensitivity data
ds = pd.read_csv(ds_fp)
ds[['GeneSymbol','HGNC_ID']] = ds['Gene Symbol /Region Name'].str.split('HGNC', expand=True)
ds['HGNC/Dosage ID'] = 'HGNC'+ds['HGNC_ID']
ds.rename(columns={'HGNC/Dosage ID':'GENE ID (HGNC)'},inplace=True)
ds.drop(labels=['HGNC_ID'], inplace=True, axis=1)

# merge with genes from D&D set
gene_df=exon_df[['hgnc_symbol']].drop_duplicates(ignore_index=True).sort_values(by='hgnc_symbol',ignore_index=True)
gd.rename(columns={'GENE SYMBOL':'hgnc_symbol'},inplace=True) # rename column for merge
ds.rename(columns={'GeneSymbol':'hgnc_symbol'},inplace=True)
gd = gd[['hgnc_symbol', 'GENE ID (HGNC)', 'DISEASE LABEL','DISEASE ID (MONDO)', 'MOI','CLASSIFICATION']] # pare down some columns
ds = ds[['hgnc_symbol', 'HI Score','%HI','pLI','LOEUF']]
im_df = gene_df.merge(gd, on='hgnc_symbol', how='left') # intermediate data frame
im_df = im_df.merge(ds, on='hgnc_symbol', how='left')

# %% [markdown] cell 9
# ## Now let's get the pathogenic variant information per gene.

# %% cell 10
# load pathogenic variant file
pv_file = CONKLIN_ROOT / "scripts/DN_GTEx/per_gene_checkpoint.csv"
pv_df=pd.read_csv(pv_file)
# merge with intermediate df
pv_df = pv_df.rename(columns={'gene':'hgnc_symbol'}).drop(labels='counter', axis=1)
im_df = im_df.merge(pv_df, on='hgnc_symbol', how='left')

# %% [markdown] cell 11
# ## Get s_het evolutionary constraint metrics for each gene

# %% cell 12
# load s_het file
s_het = pd.read_table(CONKLIN_ROOT / 'data/dHS_and_related_GeneSets/s_het/s_het_estimates.genebayes.tsv')
s_het.rename(columns={'hgnc':'GENE ID (HGNC)'},inplace=True)
im_df=im_df.merge(s_het, on='GENE ID (HGNC)', how='left')

# %% [markdown] cell 13
# ## Get pLI scores from gnomAD for each gene
# We'll do this because ClinGen often doesn't include the structured pLI information in their downloadable files.

# %% cell 14
# load gnomAD data
constraint=pd.read_table(CONKLIN_ROOT / 'data/dHS_and_related_GeneSets/constraint/gnomad.v4.1.constraint_metrics.tsv')
constraint_canonical = constraint[constraint['canonical']==True]

# convert all current pLI values to numeric
im_df['pLI'] = pd.to_numeric(im_df["pLI"])

# fill pLI field if score is missing for genes in the data frame
for gene in im_df.hgnc_symbol.unique():
    gene_filt=im_df[im_df['hgnc_symbol']==gene]
    if (math.isnan(gene_filt.pLI.values[0])):
        # get the gene's canonical transcript pLI from gnomAD
        const_filt = constraint_canonical[constraint_canonical['gene']==gene]
        if len(const_filt.gene)==0:
            const_filt = constraint_canonical[constraint_canonical['gene_id']==gene_filt.ensg.values[0]]
        if len(const_filt.gene)>0:
            cur_pLI=const_filt['lof.pLI'].values[0]
            if math.isnan(cur_pLI)==False:
                # get the location in im_df where the current gene is
                indices = im_df.index[im_df["hgnc_symbol"] == gene].tolist()
                # fill those locations with the new pLI val
                for idx in indices:
                    im_df.loc[idx,'pLI']=cur_pLI

# %% cell 15

# %% [markdown] cell 16
# ## Note the genes targetable by each editing strategy
# Note this is just based on which genes have common vars targetable by each strategy, regardless of PAMs or a CRISPR/Cas9 construct.

# %% cell 17
# set up the editing strategy names
strats=['indels','CRISPRoff','donor_base_edits','acceptor_base_edits','excision']
non_excision_fps="ubiq_region_CommonVars/CommonVars_ALL_summary.txt"
excision_fp="CommonVars/refined_common_vars/ALLgene_refined_snps.pkl"

# load the targetable gene information by strategy
all_dfs=None
for strat in strats:
    if strat=='excision':
        genes=[]
        num_vars=[]
        with open(os.path.join(results_dir, run_name, strat, excision_fp), 'rb') as fp:
            cur_dict=pickle.load(fp)
        for gene,snp_list in cur_dict.items():
            if gene!='ALLgene':
                if len(snp_list)>0:
                    if type(snp_list[0])!=type('string'):
                        genes.append(gene)
                        num_vars.append(len(snp_list))
        cur_df=pd.DataFrame({
            'gene':genes,
            'num_common_vars':num_vars,
        'strat':strat})
        cur_df_filt=cur_df[cur_df['num_common_vars']>=2]
    else:
        cur_df=pd.read_csv(os.path.join(results_dir, run_name, strat, non_excision_fps), sep='\t')
        cur_df=cur_df.iloc[:,:3]
        cur_df.columns=['index_col','gene','num_common_vars']
        cur_df.drop(labels=['index_col'],inplace=True, axis=1)
        cur_df['strat']=strat
        cur_df_filt=cur_df[cur_df['num_common_vars']>0]
    all_dfs=pd.concat([all_dfs, cur_df_filt])

# get number of genes targetable
all_genes=list(set(all_dfs.gene))

# %% cell 18

# %% cell 19
# note which gene is targetable by which strategy
# combine all of the information per strategy, to get the number of targetable genes for each strategy
indel_targetable=[] # indel targetable
indel_vars=[]
coff_targetable=[] # crisproff targetable
coff_vars=[]
be_targetable=[] # base editable at either splice donor or acceptor site
be_vars=[]
excision_targetable=[] # excision targetable
excision_vars=[]
excision_pairs=[]
gene_list=[]
for gene in list(all_dfs.gene.unique()):
    gene_list.append(gene)
    cur_gene_df=all_dfs[all_dfs.gene==gene]
    gene_strats = list(cur_gene_df.strat.values)

    # statements to store if gene is targetable by each strat
    if 'indels' in gene_strats:
        indel_targetable.append(1)
        indel_vars.append((cur_gene_df[cur_gene_df['strat']=='indels'])['num_common_vars'].values[0])
    else:
        indel_targetable.append(0)
        indel_vars.append(None)
        
    if 'CRISPRoff' in gene_strats:
        coff_targetable.append(1)
        coff_vars.append((cur_gene_df[cur_gene_df['strat']=='CRISPRoff'])['num_common_vars'].values[0])
    else:
        coff_targetable.append(0)
        coff_vars.append(None)
        
    if ('donor_base_edits' in gene_strats) or ('acceptor_base_edits' in gene_strats):
        be_targetable.append(1)
        # get number of unique vars across both of the base edit strats
        cur_vars=[]
        cur_gene_donor_vars=[]
        cur_gene_acc_vars=[]
        if ('donor_base_edits' in gene_strats):
            with open(os.path.join(results_dir, run_name, 'donor_base_edits', 'ubiq_region_CommonVars/CommonVars_ALL_dict.pkl'), 'rb') as fp:
                donor_dict=pickle.load(fp)
            cur_gene_donor_vars=donor_dict[gene]
            cur_gene_donor_vars = [x[0] for x in cur_gene_donor_vars] # separate the var position from the allele frequency
        if ('acceptor_base_edits' in gene_strats):
            with open(os.path.join(results_dir, run_name, 'acceptor_base_edits', 'ubiq_region_CommonVars/CommonVars_ALL_dict.pkl'), 'rb') as fp:
                acc_dict=pickle.load(fp)
            cur_gene_acc_vars=acc_dict[gene]
            cur_gene_acc_vars = [x[0] for x in cur_gene_acc_vars]
        be_vars_im=set(cur_gene_donor_vars+cur_gene_acc_vars)
        be_vars.append(len(be_vars_im))
    else:
        be_targetable.append(0)
        be_vars.append(None)
        
    if 'excision' in gene_strats:
        excision_targetable.append(1)
        excision_vars.append((cur_gene_df[cur_gene_df['strat']=='excision'])['num_common_vars'].values[0])
        # load valid excision pairs file
        with open(os.path.join(results_dir, run_name, 'excision/CommonVars/valid_snp_pairs/' + gene + '_valid_snp_pairs.pkl'), 'rb') as fp:
            cur_list=pickle.load(fp)
        excision_pairs.append(len(cur_list))
    else:
        excision_targetable.append(0)
        excision_vars.append(None)
        excision_pairs.append(None)

# combine into one dataframe
final_df = pd.DataFrame({
    'gene':gene_list,
    'indel_targetable':indel_targetable,
    'num_indel_vars':indel_vars,
    'crisproff_targetable':coff_targetable,
    'num_crisproff_vars':coff_vars,
    'base_editable':be_targetable,
    'num_base_editable_vars':be_vars,
    'excision_targetable':excision_targetable,
    'num_excision_vars':excision_vars,
    'num_excision_pairs':excision_pairs
})

# %% cell 20
# Note that one gene had too many excision vars to run downstream analyses efficiently on (DCC). 
# We downsampled its vars to the top 258, so let's adjust its excision numbers accordingly

# %% cell 21
# modify
idx = final_df.index[final_df['gene'] == 'DCC'][0]
final_df.loc[idx, 'num_excision_vars'] = 258

# %% cell 22
# great, now let's also adjust the number of excision pairs

# load the text file containing the downsampled snp positions
ds_snps = pd.read_table(os.path.join(results_dir, run_name, 'excision/excavate/CommonVar_locs/DCC_CommonVar_locs_ds.txt'), names=['chrom','pos'])
# load valid snp pairs for DCC
valid_pairs_fp = os.path.join(results_dir, run_name, 'excision/CommonVars/valid_snp_pairs/DCC_valid_snp_pairs.pkl')
with open(valid_pairs_fp, 'rb') as fp:
    valid_pairs=pickle.load(fp)
# see how many pairs of the ds_snps are represented in valid_pairs
valid_pairs_set = set(valid_pairs)
combos=sorted(list(combinations(ds_snps.pos, 2)))
filtered_dcc_pairs = [c for c in combos if c in valid_pairs_set]

# %% cell 23
# change the DCC entry accordingly
idx = final_df.index[final_df['gene'] == 'DCC'][0]
final_df.loc[idx, 'num_excision_pairs'] = len(filtered_dcc_pairs)

# %% cell 24
# merge this up with the existing data frame
final_df.rename(columns={'gene':'hgnc_symbol'},inplace=True)
im_df=im_df.merge(final_df, on='hgnc_symbol', how='left')

# %% [markdown] cell 25
# ## Get HPO terms and long gene names

# %% cell 26
# load table with gene HPO terms
hpo = CONKLIN_ROOT / 'data/DiseaseManifesting_Annotations/HPO/gene_top_level_terms/gene_term_df.csv'
hpo_df=pd.read_csv(hpo,index_col=0)
hpo_df.rename(columns={'gene':'hgnc_symbol'},inplace=True)
im_df = im_df.merge(hpo_df, on='hgnc_symbol',how='left')

# %% cell 27

# %% cell 28
# load long gene name information
nameInfo=pd.read_csv(DND_HGNC_DIR / "dhs_hgnc_mapped.csv")
nameInfo.rename(columns={'Input':'hgnc_symbol'},inplace=True)
nameInfo=nameInfo[nameInfo['Match type']=='Approved symbol']
im_df=im_df.merge(nameInfo[['hgnc_symbol','Approved name']], on='hgnc_symbol',how='left')

# %% cell 29
# save a checkpoint of the df here
im_df.to_csv(checkpoint_csv, index=False)

# %% [markdown] cell 30
# ## Get the number of heterozygotes targetable for each gene

# %% cell 31
# note how many heterozygotes are targetable across all editing strategies
# combine all of the information per strategy, to get the number of targetable genes for each strategy
indel_hets_all=[] # indel targetable
coff_hets_all=[] # crisproff targetable
be_hets_all=[] # base editable at either splice donor or acceptor site
dbe_hets_all=[]
abe_hets_all=[]
excision_hets_all=[] # excision targetable
all_strat_hets=[]
gene_list=[]

# load het dicts
with open(os.path.join(results_dir, run_name, 'indels', 'excavate/het_individuals/unique_hets_hit.pkl'), 'rb') as fp:
    indel_dict=pickle.load(fp)
with open(os.path.join(results_dir, run_name, 'CRISPRoff', 'excavate/het_individuals/unique_hets_hit.pkl'), 'rb') as fp:
    coff_dict=pickle.load(fp)
with open(os.path.join(results_dir, run_name, 'donor_base_edits', 'excavate/het_individuals/unique_hets_hit.pkl'), 'rb') as fp:
    dbe_dict=pickle.load(fp)
with open(os.path.join(results_dir, run_name, 'acceptor_base_edits', 'excavate/het_individuals/unique_hets_hit.pkl'), 'rb') as fp:
    abe_dict=pickle.load(fp)

for gene in list(im_df.hgnc_symbol.unique()):
    gene_list.append(gene)
    cur_gene_df=all_dfs[all_dfs.gene==gene]

    # indels
    if gene in list(indel_dict.keys()):
        indel_hets=indel_dict[gene]
    else:
        indel_hets=set()
    indel_hets_all.append(len(indel_hets))
    # crisproff
    if gene in list(coff_dict.keys()):
        coff_hets=coff_dict[gene]
    else:
        coff_hets=set()
    coff_hets_all.append(len(coff_hets))
    # donor base edits
    if gene in list(dbe_dict.keys()):
        dbe_hets=dbe_dict[gene]
    else:
        dbe_hets=set()
    dbe_hets_all.append(len(dbe_hets))
    # acceptor base edits
    if gene in list(abe_dict.keys()):
        abe_hets=abe_dict[gene]
    else:
        abe_hets=set()
    abe_hets_all.append(len(abe_hets))
    
    # # load excision hets for the gene
    # if os.path.exists(os.path.join(results_dir, run_name, 'excision', 'excavate/guide_numbers/unique_hets_hit_' + gene + '.pkl')):
    #     with open(os.path.join(results_dir, run_name, 'excision', 'excavate/guide_numbers/unique_hets_hit_' + gene + '.pkl'), 'rb') as fp:
    #         excision_hets=pickle.load(fp)
    # else:
    #     excision_hets=set()

    # new files to use
    if os.path.exists(os.path.join(results_dir, run_name, 'excision', 'het_individuals',gene+'_hets_postPAM_checkpoint.pkl')):
        with open(os.path.join(results_dir, run_name, 'excision', 'het_individuals',gene+'_hets_postPAM_checkpoint.pkl'),'rb') as fp:
            excision_hets=pickle.load(fp)
        excision_hets=excision_hets[gene]
    else:
        excision_hets=set()
    excision_hets_all.append(len(excision_hets))


    # get the number of hets across editing strategies
    unique_hets_all_strats= indel_hets | coff_hets | dbe_hets | abe_hets | excision_hets
    be_hets = dbe_hets | abe_hets
    be_hets_all.append(len(be_hets))
    all_strat_hets.append(len(unique_hets_all_strats))


num_samples=2548
# combine into one dataframe
final_df = pd.DataFrame({
    'gene':gene_list,
    'hets_across_strats':all_strat_hets,
    'prop_hets_across_strats':[x/num_samples for x in all_strat_hets],
    'num_indel_hets':indel_hets_all,
    'indel_hets_prop':[x/num_samples for x in indel_hets_all],
    'num_crisproff_hets':coff_hets_all,
    'crisproff_hets_prop':[x/num_samples for x in coff_hets_all],
    'num_base_edit_hets':be_hets_all,
    'base_edit_hets_prop':[x/num_samples for x in be_hets_all],
    'num_excision_hets':excision_hets_all,
    'excision_hets_prop':[x/num_samples for x in excision_hets_all]
})

# %% cell 34
# if print was tested, excision hets was found

# %% [markdown] cell 35
# # Add a column that denotes if the gene was targetable by CRISPR/Cas9

# %% cell 36
# add a column that denotes if the gene was targetable by CRISPR/Cas9 constructs

# load PAM-targetable genes
indel_pam_targetable=pd.read_table(os.path.join(results_dir, run_name, 'indels','excavate/het_individuals/metadata/genes_w_valid_guides.txt'), names=['gene','chrom'])
coff_pam_targetable=pd.read_table(os.path.join(results_dir, run_name, 'CRISPRoff','excavate/het_individuals/metadata/genes_w_valid_guides.txt'),names=['gene','chrom'])
dbe_pam_targetable=pd.read_table(os.path.join(results_dir, run_name, 'donor_base_edits','excavate/het_individuals/metadata/genes_w_valid_guides.txt'),names=['gene','chrom'])
abe_pam_targetable=pd.read_table(os.path.join(results_dir, run_name, 'acceptor_base_edits','excavate/het_individuals/metadata/genes_w_valid_guides.txt'),names=['gene','chrom'])
excision_pam_targetable=pd.read_table(os.path.join(results_dir, run_name, 'excision','excavate/het_individuals/metadata/genes_w_valid_guides.txt'),names=['gene','chrom'])

# note each gene in the df
ipt=[]
cpt=[]
dpt=[]
apt=[]
ept=[]
gene_list=[]
for gene in im_df.hgnc_symbol.unique():
    gene_list.append(gene)
    if gene in list(indel_pam_targetable.gene):
        ipt.append(1)
    else:
        ipt.append(0)
    if gene in list(coff_pam_targetable.gene):
        cpt.append(1)
    else:
        cpt.append(0)
    if gene in list(dbe_pam_targetable.gene):
        dpt.append(1)
    else:
        dpt.append(0)
    if gene in list(abe_pam_targetable.gene):
        apt.append(1)
    else:
        apt.append(0)
    if gene in list(excision_pam_targetable.gene):
        ept.append(1)
    else:
        ept.append(0)
im_be2=[]
im_be=[x + y for x, y in zip(dpt, apt)]
for x in im_be:
    if x>0:
        im_be2.append(1)
    else:
        im_be2.append(0)
        
final_df2 = pd.DataFrame({
    'gene':gene_list,
    'indel_pam_targetable':ipt,
    'crisproff_pam_targetable':cpt,
    'base_edit_pam_targetable':im_be2,
    'excision_pam_targetable':ept
})

# %% cell 37
# merge with intermediate df
final_df.rename(columns={'gene':'hgnc_symbol'},inplace=True)
final_df2.rename(columns={'gene':'hgnc_symbol'},inplace=True)
im_df=im_df.merge(final_df2,on='hgnc_symbol', how='left')
im_df=im_df.merge(final_df,on='hgnc_symbol',how='left')

# %% cell 38
# save checkpoint
im_df.to_csv(checkpoint_csv, index=False)


# %% cell 44
indel_hets_all=[] # indel targetable
coff_hets_all=[] # crisproff targetable
be_hets_all=[] # base editable at either splice donor or acceptor site
dbe_hets_all=[]
abe_hets_all=[]
excision_hets_all=[] # excision targetable
all_strat_hets=[]
gene_list=[]

# load het dicts
with open(os.path.join(results_dir, run_name, 'indels', 'prePAM_hets/unique_hets_hit_prePAM_filter.pkl'), 'rb') as fp:
    indel_dict=pickle.load(fp)
with open(os.path.join(results_dir, run_name, 'CRISPRoff', 'prePAM_hets/unique_hets_hit_prePAM_filter.pkl'), 'rb') as fp:
    coff_dict=pickle.load(fp)
with open(os.path.join(results_dir, run_name, 'donor_base_edits', 'prePAM_hets/unique_hets_hit_prePAM_filter.pkl'), 'rb') as fp:
    dbe_dict=pickle.load(fp)
with open(os.path.join(results_dir, run_name, 'acceptor_base_edits', 'prePAM_hets/unique_hets_hit_prePAM_filter.pkl'), 'rb') as fp:
    abe_dict=pickle.load(fp)

for gene in list(im_df.hgnc_symbol.unique()):
    gene_list.append(gene)
    cur_gene_df=all_dfs[all_dfs.gene==gene]

    # indels
    if gene in list(indel_dict.keys()):
        indel_hets=indel_dict[gene]
    else:
        indel_hets=set()
    indel_hets_all.append(len(indel_hets))
    # crisproff
    if gene in list(coff_dict.keys()):
        coff_hets=coff_dict[gene]
    else:
        coff_hets=set()
    coff_hets_all.append(len(coff_hets))
    # donor base edits
    if gene in list(dbe_dict.keys()):
        dbe_hets=dbe_dict[gene]
    else:
        dbe_hets=set()
    dbe_hets_all.append(len(dbe_hets))
    # acceptor base edits
    if gene in list(abe_dict.keys()):
        abe_hets=abe_dict[gene]
    else:
        abe_hets=set()
    abe_hets_all.append(len(abe_hets))
    
    # load excision hets for the gene
    excision_path=os.path.join(results_dir, run_name, 'excision', 'prePAM_hets', gene + '_hets_prePAM_checkpoint.pkl')
    if os.path.exists(excision_path):
        with open(excision_path, 'rb') as fp:
            excision_hets=pickle.load(fp)
        excision_hets=excision_hets[gene]
    else:
        excision_hets=set()
    excision_hets_all.append(len(excision_hets))


    # get the number of hets across editing strategies
    unique_hets_all_strats= indel_hets | coff_hets | dbe_hets | abe_hets | excision_hets
    be_hets = dbe_hets | abe_hets
    be_hets_all.append(len(be_hets))
    all_strat_hets.append(len(unique_hets_all_strats))


num_samples=2548
# combine into one dataframe
final_df = pd.DataFrame({
    'gene':gene_list,
    'hets_across_strats_prePAM':all_strat_hets,
    'prop_hets_across_strats_prePAM':[x/num_samples for x in all_strat_hets],
    'num_indel_hets_prePAM':indel_hets_all,
    'indel_hets_prop_prePAM':[x/num_samples for x in indel_hets_all],
    'num_crisproff_hets_prePAM':coff_hets_all,
    'crisproff_hets_prop_prePAM':[x/num_samples for x in coff_hets_all],
    'num_base_edit_hets_prePAM':be_hets_all,
    'base_edit_hets_prop_prePAM':[x/num_samples for x in be_hets_all],
    'num_excision_hets_prePAM':excision_hets_all,
    'excision_hets_prop_prePAM':[x/num_samples for x in excision_hets_all]
})

# %% cell 45
# combine with existing df
final_df.rename(columns={'gene':'hgnc_symbol'},inplace=True)
im_df=im_df.merge(final_df,on='hgnc_symbol', how='left')

# %% cell 46
# save checkpoint
im_df.to_csv(checkpoint_csv, index=False)

# %% [markdown] cell 47
# ## Reload the data frame here

# %% cell 48
im_df = pd.read_csv(checkpoint_csv)

# %% [markdown] cell 49
# ## Get guide information

# %% [markdown] cell 50
# ### Number of guides required to target all possible heterozygotes

# %% cell 51
# load guide information for each editing strategy (except excision, we'll do that separately)
ig = pd.read_csv(os.path.join(results_dir, run_name, 'indels', 'excavate/guide_numbers/guide_info_by_gene.csv'),index_col=0)
ig['strat']='indels'
cg = pd.read_csv(os.path.join(results_dir, run_name, 'CRISPRoff', 'excavate/guide_numbers/guide_info_by_gene.csv'),index_col=0)
cg['strat']='CRISPRoff'
bg = pd.read_csv(os.path.join(results_dir, run_name, 'summary_files/combined_base_edit_guide_info.csv'), index_col=0)
bg['strat']='base_edits'
all_dfs = pd.concat([ig,cg,bg])

# %% cell 52
excision_guide_fp="excision/excavate/guide_numbers"
filenames=os.listdir(os.path.join(results_dir, run_name, excision_guide_fp))
filenames=[x for x in filenames if 'guide_info_by_gene.csv' not in x]
filenames=[x for x in filenames if 'no_hets_notes' not in x]
filenames=[x for x in filenames if 'unique_hets_hit' not in x]
filenames.sort()
excision_dfs=None
for file in filenames:
    excision_dfs=pd.concat([excision_dfs,pd.read_csv(os.path.join(results_dir, run_name, excision_guide_fp, file))])
excision_dfs['strat']='excision'
all_dfs = pd.concat([all_dfs,excision_dfs])

# %% cell 53
strats=['indels','CRISPRoff','base_edits', 'excision']
from ast import literal_eval
gene_dfs=None
for strat in strats:
    strat_df=all_dfs[all_dfs['strat']==strat]
    for gene in strat_df.gene.unique():
        gene_df=strat_df[strat_df['gene']==gene]
        guide_set=set()
        num_unique_guides=[]
        if strat=='excision':
            for pair in gene_df.selected_snps:
                pair=literal_eval(pair)
                if (pair[0] in guide_set)==False:
                    guide_set.add(pair[0])
                if (pair[1] in guide_set)==False:
                    guide_set.add(pair[1])
                num_unique_guides.append(len(guide_set))
        else:
            for snp in gene_df.selected_snps:
                if (snp in guide_set)==False:
                    guide_set.add(snp)
                num_unique_guides.append(len(guide_set))
        gene_df['num_unique_guides']=num_unique_guides
        gene_dfs=pd.concat([gene_dfs,gene_df])

# %% cell 54

# %% cell 55
# get number of hets targeted by the first four unique gRNAs
four_guides = gene_dfs[gene_dfs['num_unique_guides']<=4]
num_hets = pd.DataFrame(four_guides.groupby(['gene','strat'])['selected_hets_steps'].max())
num_hets.reset_index(inplace=True)
num_hets.drop_duplicates(inplace=True, ignore_index=True)

# %% cell 56
# append this to the data frame
ih =[]
ig=[]
ch=[]
cg=[]
bh=[]
bg=[]
eh=[]
eg=[]
for idx,row in num_hets.iterrows():
    if row.strat=='indels':
        ig.append(row.gene)
        ih.append(row.selected_hets_steps)
    if row.strat=='CRISPRoff':
        cg.append(row.gene)
        ch.append(row.selected_hets_steps)
    if row.strat=='base_edits':
        bg.append(row.gene)
        bh.append(row.selected_hets_steps)
    if row.strat=='excision':
        eg.append(row.gene)
        eh.append(row.selected_hets_steps)
idf=pd.DataFrame({'gene':ig,'num_hets_targeted_four_guides_indels':ih})
cdf=pd.DataFrame({'gene':cg,'num_hets_targeted_four_guides_crisproff':ch})
bdf=pd.DataFrame({'gene':bg,'num_hets_targeted_four_guides_base_edits':bh})
edf=pd.DataFrame({'gene':eg,'num_hets_targeted_four_guides_excision':eh})
combined_df=idf.merge(cdf, on='gene', how='outer')
combined_df=combined_df.merge(bdf, on='gene',how='outer')
combined_df=combined_df.merge(edf, on='gene',how='outer')
num_samples=2548
combined_df['prop_hets_targeted_four_guides_indels']=combined_df['num_hets_targeted_four_guides_indels']/2548
combined_df['prop_hets_targeted_four_guides_crisproff']=combined_df['num_hets_targeted_four_guides_crisproff']/2548
combined_df['prop_hets_targeted_four_guides_base_edits']=combined_df['num_hets_targeted_four_guides_base_edits']/2548
combined_df['prop_hets_targeted_four_guides_excision']=combined_df['num_hets_targeted_four_guides_excision']/2548

# %% cell 57
# merge up with larger df
combined_df.rename(columns={'gene':'hgnc_symbol'}, inplace=True)
im_df = im_df.merge(combined_df, on='hgnc_symbol', how='left')

# %% cell 59
# save checkpoint
im_df.to_csv(checkpoint_csv, index=False)

# %% [markdown] cell 60
# # Get the number of genes that are indel targetable post-NMD filtering

# %% cell 61
nmd_summary=pd.read_csv(os.path.join(results_dir, run_name, "indels/NMD/NMD_induction_summary.csv"))
nmd_vars = pd.read_csv(os.path.join(results_dir, run_name, "indels/NMD/NMD_induction_var_info.csv"))

# %% cell 64
num_vars_inducing_nmd = []
gene_list = []
indel_targetable_post_NMD_assessment=[]
for idx,row in nmd_vars.iterrows():
    gene_list.append(row.gene)
    vars_inducing_NMD=literal_eval(row.vars_consistently_inducing_NMD)
    if len(vars_inducing_NMD)>0:
        num_vars_inducing_nmd.append(len(vars_inducing_NMD))
        indel_targetable_post_NMD_assessment.append(1)
    else:
        indel_targetable_post_NMD_assessment.append(0)
        num_vars_inducing_nmd.append(0)
nmd_df = pd.DataFrame({
    'hgnc_symbol':gene_list,
    'indel_targetable_post_NMD_assessment':indel_targetable_post_NMD_assessment,
    'num_vars_inducing_NMD':num_vars_inducing_nmd
})

# %% cell 65
# merge back with intermediate data frame
im_df = im_df.merge(nmd_df, on='hgnc_symbol', how='left')

# %% cell 66
im_df = im_df[im_df['MOI'].isin(['AD','SD'])]

# %% cell 67
# rename some columns
im_df.rename(columns={'indel_targetable':'indel_targetable_pre_NMD_assessment'},inplace=True)

# %% cell 68
# save a checkpoint
im_df.to_csv(checkpoint_csv, index=False)

# %% [markdown] cell 70
# # Get the number of guides a different way-multiple strategy prioritization

# %% cell 71
# Reuse the selected run configuration when reading cross-strategy guide summaries.
summary_dir="summary_files/cross_strat_gRNAs/excision/results"

files = os.listdir(os.path.join(results_dir, run_name, summary_dir))
full_df=None
for file in files:
    cur_gene = file.split('_')[0]
    cur_df=pd.read_csv(os.path.join(results_dir, run_name, summary_dir, file),index_col=0)
    cur_df['gene']=cur_gene
    full_df = pd.concat([full_df,cur_df])

guide_cutoff=4
gene_list=[]
num_hets_four_guides=[]
num_hets_all_guides=[]
num_haps_four_guides=[]
num_haps_all_guides=[]
max_guides=[]
for gene in full_df.gene.unique():
    gene_list.append(gene)
    filt_df = full_df[full_df['gene']==gene]
    four_df = filt_df[filt_df['num_guides']<=guide_cutoff]

    # update numbers for four-guide limit
    num_hets_four_guides.append(max(four_df.cumulative_people_targeted))
    num_haps_four_guides.append(max(four_df.cumulative_haplotypes_targeted))

    # update numbers for no-guide limit (with max iterations at 50)
    num_hets_all_guides.append(max(filt_df.cumulative_people_targeted))
    num_haps_all_guides.append(max(filt_df.cumulative_haplotypes_targeted))

    # update guide numbers
    max_guides.append(max(filt_df.num_guides))
    
excision_guide_info=pd.DataFrame({
    'gene':gene_list,
    'excision_num_haps_four_guides_alg2':num_haps_four_guides,
    'excision_num_hets_four_guides_alg2':num_hets_four_guides,
    'excision_num_haps_all_guides_alg2':num_haps_all_guides,
    'excision_num_hets_all_guides_alg2':num_hets_all_guides,
    'excision_max_guides':max_guides
})

# %% cell 72
# non-excision guides
summary_dir = "summary_files/cross_strat_gRNAs/non_excision_strats_separate/results"
strats=['indels','CRISPRoff','base_edit_', 'donor_base_edits','acceptor_base_edits'] # need to switch this to just 'base edits'

# get the per-strategy people & haplotypes targeted for the non-excision strats
files=os.listdir(os.path.join(results_dir, run_name, summary_dir))
by_strat_df=None
for strat in strats:
    cur_strat_files=[x for x in files if strat in x]
    for file in cur_strat_files:
        cur_gene = file.split('_')[0]
        cur_df = pd.read_csv(os.path.join(results_dir, run_name, summary_dir,file),index_col=0)
        cur_df['strat']=strat
        cur_df['gene']=cur_gene
        by_strat_df=pd.concat([by_strat_df,cur_df])

# same as before, see how many haps and people are captured under various filters
guide_cutoff=4
strat_dict={}
for strat in strats:
    gene_list=[]
    num_hets_four_guides=[]
    num_hets_all_guides=[]
    num_haps_four_guides=[]
    num_haps_all_guides=[]
    max_guides=[]
    
    for gene in by_strat_df.gene.unique():
        filt_df = by_strat_df[by_strat_df['gene']==gene]
        filt_df = filt_df[filt_df['strat']==strat]
        if len(filt_df)==0:
            continue
        gene_list.append(gene)
        four_df = filt_df[filt_df['num_guides']<=guide_cutoff]
            
    
        # update numbers for four-guide limit
        num_hets_four_guides.append(max(four_df.cumulative_people_targeted))
        num_haps_four_guides.append(max(four_df.cumulative_haplotypes_targeted))
    
        # update numbers for no-guide limit (with max iterations at 50)
        num_hets_all_guides.append(max(filt_df.cumulative_people_targeted))
        num_haps_all_guides.append(max(filt_df.cumulative_haplotypes_targeted))
        max_guides.append(max(filt_df.num_guides))
        
    ne_guide_info=None
    ne_guide_info=pd.DataFrame({
        'gene':gene_list,
        strat + '_num_haps_four_guides_alg2':num_haps_four_guides,
        strat + '_num_hets_four_guides_alg2':num_hets_four_guides,
        strat + '_num_haps_all_guides_alg2':num_haps_all_guides,
        strat + '_num_hets_all_guides_alg2':num_hets_all_guides,
        strat + '_max_guides': max_guides
    })
    strat_dict[strat]=ne_guide_info

# %% cell 73
# merge these all with the final data frame
for strat,cur_df in strat_dict.items():
    cur_df.rename(columns={'gene':'hgnc_symbol'},inplace=True)
    im_df = im_df.merge(cur_df, on='hgnc_symbol',how='left')
excision_guide_info.rename(columns={'gene':'hgnc_symbol'},inplace=True)
im_df = im_df.merge(excision_guide_info, how='left', on='hgnc_symbol')

# %% cell 74
# note that there will be repeat information in the acceptor, donor, and general base edit columns

# %% cell 75

# %% cell 76
metrics = [
    c.replace('base_edit__', '')
    for c in im_df.columns
    if c.startswith('base_edit__')
]
for m in metrics:
    base_col = f'base_edit__{m}'
    donor_col = f'donor_base_edits_{m}'
    acceptor_col = f'acceptor_base_edits_{m}'
    
    im_df[f'combined_base_edit_{m}'] = (
        im_df[base_col]
            .combine_first(im_df[donor_col])
            .combine_first(im_df[acceptor_col])
    )

# %% cell 77
# clean this up a bit
im_df.drop(labels=['base_edit__num_haps_four_guides_alg2',
       'base_edit__num_hets_four_guides_alg2',
       'base_edit__num_haps_all_guides_alg2',
       'base_edit__num_hets_all_guides_alg2', 'base_edit__max_guides',
       'donor_base_edits_num_haps_four_guides_alg2',
       'donor_base_edits_num_hets_four_guides_alg2',
       'donor_base_edits_num_haps_all_guides_alg2',
       'donor_base_edits_num_hets_all_guides_alg2',
       'donor_base_edits_max_guides',
       'acceptor_base_edits_num_haps_four_guides_alg2',
       'acceptor_base_edits_num_hets_four_guides_alg2',
       'acceptor_base_edits_num_haps_all_guides_alg2',
       'acceptor_base_edits_num_hets_all_guides_alg2',
       'acceptor_base_edits_max_guides'],
           axis=1,inplace=True)

# %% cell 78
# clean up further
im_df.drop(labels=['num_hets_targeted_four_guides_crisproff',
       'num_hets_targeted_four_guides_base_edits',
       'num_hets_targeted_four_guides_excision',
       'prop_hets_targeted_four_guides_indels',
       'prop_hets_targeted_four_guides_crisproff',
       'prop_hets_targeted_four_guides_base_edits',
       'prop_hets_targeted_four_guides_excision'],axis=1,inplace=True)

# %% cell 79
# great, now save this
im_df.to_csv(os.path.join(results_dir, run_name, 'summary_files/dnd_593.csv'),index=False)

# %% [markdown] cell 81
# ## Now, since we want to consider each set of excision guides as 1 'therapy', let's condense some of the numbers

# %% cell 82
# reload

# %% cell 83
# Reuse the selected run configuration when reloading cross-strategy guide summaries.
summary_dir="summary_files/cross_strat_gRNAs/excision/results"

files = os.listdir(os.path.join(results_dir, run_name, summary_dir))
im_df = pd.read_csv(checkpoint_csv)
full_df=None
for file in files:
    cur_gene = file.split('_')[0]
    cur_df=pd.read_csv(os.path.join(results_dir, run_name, summary_dir, file),index_col=0)
    cur_df['gene']=cur_gene
    full_df = pd.concat([full_df,cur_df])

# %% cell 85
# create a 'therapies' column
full_df['therapies'] = (
    full_df.groupby('gene')['num_guides']
      .diff()                      # change in guides within gene
      .gt(0)                       # True when guides increase
      .groupby(full_df['gene'])
      .cumsum()                    # cumulative therapy count
      .add(1)                      # start therapy numbering at 1
)

# %% cell 86
# get the cutoff at 4 therapies for excisions
therapy_cutoff=4
gene_list=[]
num_hets_four_guides=[]
num_haps_four_guides=[]
max_therapies=[]
for gene in full_df.gene.unique():
    gene_list.append(gene)
    filt_df = full_df[full_df['gene']==gene]
    four_df = filt_df[filt_df['therapies']<=therapy_cutoff]

    # update numbers for four-guide limit
    num_hets_four_guides.append(max(four_df.cumulative_people_targeted))
    num_haps_four_guides.append(max(four_df.cumulative_haplotypes_targeted))

    # update guide numbers
    max_therapies.append(max(filt_df.num_guides))
    
excision_guide_info=pd.DataFrame({
    'gene':gene_list,
    'excision_num_haps_four_therapies_alg2':num_haps_four_guides,
    'excision_num_hets_four_therapies_alg2':num_hets_four_guides,
    'excision_max_therapies':max_therapies
})

# %% cell 88
excision_guide_info.rename(columns={'gene':'hgnc_symbol'},inplace=True)
im_df = im_df.merge(excision_guide_info, on='hgnc_symbol',how='left')

# %% cell 89
# save the df
Path(output_csv).parent.mkdir(parents=True, exist_ok=True)
im_df.to_csv(output_csv, index=False)
