# imports
import os
import sys
import numpy as np
import pandas as pd
from datetime import datetime
from itertools import permutations, combinations, product
import itertools
#import seaborn as sns
from matplotlib.ticker import FuncFormatter
import pickle
import math
from ast import literal_eval
import gzip
import argparse
import time
import pysam

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--results_dir', type = str, required = True, help = 'Directory with the results of editing strategies and common variants.')
    parser.add_argument('--run_name',type=str,required=True,help='run name for the results directory.')
    parser.add_argument('--gene', type=str,required=True,help="Current gene name.")
    parser.add_argument('--ref_genome_fasta', type=str, required=True, help='Reference genome fasta file path.')
    parser.add_argument('--bt_dir', type=str, required=True, help='Browswer track results directory')
    parser.add_argument('--exon_file', type=str, required=True, help='Exon information across genes.')
    parser.add_argument('--sample_map', type=str, required=True, help='Filepath mapping sampled IDs to more specific sample names.')
    parser.add_argument('--rsID_fp', type=str, required=True, help='map between position and rsID')
    args = parser.parse_args()
    return args


def main():
    # parse input args
    args=parse_args()
    results_dir=args.results_dir
    run_name=args.run_name
    gene=args.gene
    ref_genome_fasta=args.ref_genome_fasta
    exon_file=args.exon_file
    sample_map_fp=args.sample_map
    rsID_fp=args.rsID_fp
    strats=['indels','CRISPRoff','donor_base_edits', 'acceptor_base_edits','excision']


    ## Load pre-PAM SNPs for each gene
    def normalize_pos_af_list(raw_vals):
        if not isinstance(raw_vals, list) or len(raw_vals) == 0:
            return []
        if not isinstance(raw_vals[0], list):
            return []
        return raw_vals

    def coerce_pos_numeric(df, pos_col='pos'):
        if pos_col not in df.columns:
            return df
        df = df.copy()
        df[pos_col] = pd.to_numeric(df[pos_col], errors='coerce')
        return df[df[pos_col].notna()].copy()

    pre_pam_frames = []
    for strat in strats:
        # load any excision snps
        if strat=='excision':
            excision_df = None
            excision_path = os.path.join(results_dir, run_name, strat, 'CommonVars', 'refined_common_vars',gene + '_refined_snp_list.pkl')
            if os.path.exists(excision_path):
                with open(excision_path, 'rb') as fp:
                    excision_snps = pickle.load(fp)
                if len(excision_snps)>0:
                    excision_df = pd.DataFrame(sorted(excision_snps), columns=['pos'])
                    # load the file with allele frequencies
                    with open(os.path.join(results_dir, run_name, strat, 'CommonVars','CommonVars_ALL_dict.pkl'), 'rb') as fp:
                        af_info = pickle.load(fp)
                    # extract allele info for the relevant genes
                    cur_gene_afs = af_info.get(gene, [])
                    # filter to just the positions in excision_snps set
                    filtered_pos = [x for x in cur_gene_afs if x[0] in excision_snps]
                    # now convert to df
                    excision_df = pd.DataFrame(filtered_pos, columns=['pos','af'])
                    excision_df['edit_strat'] = strat
                    pre_pam_frames.append(excision_df)
                
        # load any indel snps respecting NMD induction
        elif strat=='indels':
            indels_df = None
            indels_path = os.path.join(results_dir, run_name, strat, 'NMD', 'NMD_induction_var_info.csv')
            if os.path.exists(indels_path):
                nmd_df = pd.read_csv(indels_path)
                if gene in list(nmd_df['gene']):
                    cur_gene_nmd = nmd_df[nmd_df['gene']==gene]
                    cur_gene_nmd_vars = literal_eval(cur_gene_nmd.vars_consistently_inducing_NMD.values[0])
                    if len(cur_gene_nmd_vars)>0:
                        af_included=os.path.join(results_dir, run_name, strat, 'ubiq_region_CommonVars', 'CommonVars_ALL_dict.pkl')
                        with open(af_included, 'rb') as fp:
                            afs = pickle.load(fp)
                        cur_gene_afs = normalize_pos_af_list(afs.get(gene, []))
                        if len(cur_gene_afs)>0:
                            # filter snp list to only those that induce nmd
                            filtered_snps = [sub for sub in cur_gene_afs if sub[0] in set(cur_gene_nmd_vars)]
                            if len(filtered_snps)>0:
                                indel_df = pd.DataFrame(filtered_snps, columns=['pos','af'])
                                indel_df['edit_strat'] = strat
                                pre_pam_frames.append(indel_df)

                        
        # load any non-indel, non-excision snps 
        else:
            non_excision_df = None
            non_excision_path = os.path.join(results_dir, run_name, strat, 'ubiq_region_CommonVars', 'CommonVars_ALL_dict.pkl')
            if os.path.exists(non_excision_path):
                with open(non_excision_path, 'rb') as fp:
                    non_excision_snps = pickle.load(fp)
                cur_gene_snps = normalize_pos_af_list(non_excision_snps.get(gene, []))
                if len(cur_gene_snps)>0:
                    non_excision_df = pd.DataFrame(cur_gene_snps, columns=['pos','af'])
                    non_excision_df['edit_strat']=strat
                    pre_pam_frames.append(non_excision_df)
    if len(pre_pam_frames) > 0:
        final_df = pd.concat(pre_pam_frames, ignore_index=True)
    else:
        final_df = pd.DataFrame(columns=['pos','af','edit_strat'])
    final_df = coerce_pos_numeric(final_df, 'pos')
    # condense editing strategies
    if final_df.empty:
        prePAM_snps = pd.DataFrame(columns=['pos','af','edit_strat'])
    else:
        im = pd.DataFrame(final_df.groupby('pos')['edit_strat'].agg(set)).reset_index()
        prePAM_snps = im.merge(final_df[['pos','af']], on='pos', how='left')
        prePAM_snps['edit_strat'] = (
            prePAM_snps['edit_strat']
            .apply(lambda x: "{" + ", ".join(f"'{v}'" for v in sorted(set(x))) + "}")
        )
        prePAM_snps.drop_duplicates(inplace=True)

    # ensure that these are all biallelic--if not, remove them from consideration

    # vcf formatting info
    num_samples=2548
    sample_list = list(range(1,num_samples+1)) # number of people in 1KG
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]

    def assert_unique_and_biallelic_vcf_values(vcf):
        # 1️⃣ remove fully duplicated rows
        vcf = vcf.drop_duplicates()
        
        # 2️⃣ keep only rows where ref and alt are single letters
        vcf = vcf[vcf['ref'].str.len() == 1]
        vcf = vcf[vcf['alt'].str.len() == 1]
        
        # optional: reset index
        vcf = vcf.reset_index(drop=True)
        
        return vcf

    passing_snps=[]
    snp_pos=[]
    alt=[]
    ref=[]
    for strat in strats:
        
        vcf_fp = os.path.join(results_dir, run_name, strat, "excavate/input_vcfs/" + gene + '_CommonVar_filtered.vcf.gz')
        if os.path.exists(vcf_fp) is False:
            continue
        vcf = pd.read_table(vcf_fp, comment='#', header=None, compression='gzip')
        vcf.columns=cols
        vcf = assert_unique_and_biallelic_vcf_values(vcf)
        passing_snps.extend(vcf.pos)
        alt.extend(vcf.alt)
        ref.extend(vcf.ref)

    prePAM_snps = prePAM_snps[prePAM_snps['pos'].isin(passing_snps)]
    ref_alt_df = pd.DataFrame({'pos':passing_snps,'ref':ref,'alt':alt})
    ref_alt_df = coerce_pos_numeric(ref_alt_df, 'pos')
    prePAM_snps = prePAM_snps.merge(ref_alt_df, on='pos',how='left')
    prePAM_snps = coerce_pos_numeric(prePAM_snps, 'pos')

    ## Load post-PAM SNPS
    post_pam_frames = []
    # vcf formatting info
    num_samples=2548
    sample_list = list(range(1,num_samples+1)) # number of people in 1KG
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]

    def assert_unique_and_biallelic_vcf_values(vcf):
        # 1️⃣ remove fully duplicated rows
        vcf = vcf.drop_duplicates()
        
        # 2️⃣ keep only rows where ref and alt are single letters
        vcf = vcf[vcf['ref'].str.len() == 1]
        vcf = vcf[vcf['alt'].str.len() == 1]
        
        # optional: reset index
        vcf = vcf.reset_index(drop=True)
        
        return vcf

    for strat in strats:
        
        vcf_fp = os.path.join(results_dir, run_name, strat, "excavate/Guide_filtered_vcfs/" + gene + '_guide_filtered.vcf')
        if os.path.exists(vcf_fp) is False:
            continue
        vcf = pd.read_table(vcf_fp, comment='#', header=None)
        vcf.columns=cols
        vcf = assert_unique_and_biallelic_vcf_values(vcf)
        passing = prePAM_snps[prePAM_snps['pos'].isin(list(vcf.pos))]
        if not passing.empty:
            post_pam_frames.append(passing)

    if len(post_pam_frames) > 0:
        final_df = pd.concat(post_pam_frames, ignore_index=True)
    else:
        final_df = pd.DataFrame(columns=prePAM_snps.columns)

    if final_df.empty:
        postPAM_snps = pd.DataFrame(columns=prePAM_snps.columns)
    else:
        postPAM_snps = final_df.drop_duplicates().reset_index(drop=True)
    postPAM_snps = coerce_pos_numeric(postPAM_snps, 'pos')

    # get genomic coordinates of first and last pre-PAM SNP
    window_buffer = 500
    window_start = prePAM_snps.pos.min() - window_buffer
    window_end = prePAM_snps.pos.max() + window_buffer

    # get genotype frequencies
    # vcf formatting info
    num_samples=2548
    sample_list = list(range(1,num_samples+1)) # number of people in 1KG
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]

    def get_het_freq(cur_snp_inds):
        # remove non-sample columns
        cur_snp_inds = cur_snp_inds.iloc[:, 9:].copy()
        return het_freq, het_num

    # load vcf for each strategy
    vcf_combined = []
    for strat in strats:
        vcf_fp = os.path.join(results_dir, run_name, strat, "excavate/input_vcfs/" + gene + '_CommonVar_filtered.vcf.gz')
        if os.path.exists(vcf_fp) is False:
            continue
        cur_vcf = pd.read_table(vcf_fp, comment='#', header=None, compression='gzip')
        cur_vcf.columns=cols
        vcf_combined.append(cur_vcf)
    if len(vcf_combined) > 0:
        vcf_combined = pd.concat(vcf_combined, ignore_index=True)
        vcf_combined.drop_duplicates(ignore_index=True, inplace=True)
        vcf_combined = assert_unique_and_biallelic_vcf_values(vcf_combined)
        vcf_filt = vcf_combined[vcf_combined.pos.isin(list(prePAM_snps.pos))]
    else:
        vcf_filt = pd.DataFrame(columns=cols)

    het_freq=[]
    het_num=[]
    homo_ref=[]
    homo_ref_num=[]
    homo_alt_freq=[]
    homo_alt_num=[]
    snps=[]

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

    def make_list_of_ref_homs(r_clean_tdf):
        # for each SNP, make a list of heterozygous individuals
        indviduals = r_clean_tdf.index.tolist()
        bool_df = r_clean_tdf.map(lambda x: ((x == "0|0"))).copy()

        snp_het_indviduals_dict = dict() # make a dictionary
        for snp in r_clean_tdf.columns:
            het_individual = bool_df.index[bool_df[snp]].tolist()
            snp_het_indviduals_dict[snp] = het_individual

        return snp_het_indviduals_dict
        
    def make_list_of_alt_homs(r_clean_tdf):
        # for each SNP, make a list of heterozygous individuals
        indviduals = r_clean_tdf.index.tolist()
        bool_df = r_clean_tdf.map(lambda x: ((x == "1|1"))).copy()

        snp_het_indviduals_dict = dict() # make a dictionary
        for snp in r_clean_tdf.columns:
            het_individual = bool_df.index[bool_df[snp]].tolist()
            snp_het_indviduals_dict[snp] = het_individual

        return snp_het_indviduals_dict

    if vcf_filt.empty:
        het_dict = {}
        ref_hom_dict = {}
        alt_hom_dict = {}
    else:
        # transpose and clean the vcf
        vcf_t, key_df = transpose_and_clean(vcf_filt)
        # calculate genotype frequencies
        het_dict=make_list_of_hets(vcf_t)
        ref_hom_dict=make_list_of_ref_homs(vcf_t)
        alt_hom_dict=make_list_of_alt_homs(vcf_t)
    # set up variables for storage
    snp_list=[]
    het_freq=[]
    het_num=[]
    homo_ref_freq=[]
    homo_ref_num=[]
    homo_alt_freq=[]
    homo_alt_num=[]
    # combine the data into general frequencies
    for snp,het_list in het_dict.items():
        snp_list.append(snp)
        het_num.append(len(het_list))
        het_freq.append(len(het_list)/num_samples)
        homo_ref_num.append(len(ref_hom_dict[snp]))
        homo_ref_freq.append(len(ref_hom_dict[snp])/num_samples)
        homo_alt_num.append(len(alt_hom_dict[snp]))
        homo_alt_freq.append(len(alt_hom_dict[snp])/num_samples)

    final_df = pd.DataFrame({
        'pos':snp_list,
        'heterozygote_freq':het_freq,
        'het_num':het_num,
        'homo_ref_freq':homo_ref_freq,
        'homo_ref_num':homo_ref_num,
        'homo_alt_freq':homo_alt_freq,
        'homo_alt_num':homo_alt_num
    })
    # merge with the pre- and post-PAM dfs
    prePAM_snps = prePAM_snps.merge(final_df, on='pos', how='left')
    postPAM_snps = postPAM_snps.merge(final_df, on='pos', how='left')

    ## Get list of valid partner SNPs for excision snps
    valid_pairs_fp = os.path.join(results_dir, run_name, "excision/CommonVars/valid_snp_pairs",  gene + "_valid_snp_pairs.pkl")
    snp_partners={}
    closest_partners={}
    if os.path.exists(valid_pairs_fp):
        # load valid excision snp pairs
        with open(valid_pairs_fp,'rb') as fp:
            valid_pairs = pickle.load(fp)
        # loop through flattened SNPs and get each one's partners from the valid pairs list
        for idx,row in prePAM_snps.iterrows():
            snp=row.pos
            # only if the snp is an excision snp should we be assessing this
            if 'excision' in row.edit_strat:
                pairs_containing_snp = {x for x in valid_pairs if (x[0]==snp or x[1]==snp)}
                flat_set = {x for tup in pairs_containing_snp for x in tup}
                # remove the snp of interest from the 'partners' set
                flat_set.discard(snp)
                # store, and we can add to data frame at a later time
                snp_partners[int(snp)] = list(flat_set)
                # also store the 5 closest excision partners, as this may be easier to visualize
                ordered_partners=list(flat_set.copy())
                ordered_partners.sort()
                if len(ordered_partners)<5:
                    closest_partners[int(snp)]=ordered_partners[0:len(ordered_partners)]
                else:
                    closest_partners[int(snp)]=ordered_partners[0:5]

    ## Get gRNA sequences
    # combine excavate output files across strats
    gRNA_df = []
    for strat in strats:
        guides_path = os.path.join(results_dir, run_name, strat, "excavate", "excavate_outputs", gene+"_output", "all_guides.csv")
        if os.path.exists(guides_path):
            guide_info = pd.read_csv(guides_path)
            gRNA_df.append(guide_info)
    if len(gRNA_df) > 0:
        gRNA_df = [x for x in gRNA_df if not x.empty]
        if len(gRNA_df) > 0:
            gRNA_df = pd.concat(gRNA_df, ignore_index=True).drop_duplicates(ignore_index=True)
        else:
            gRNA_df = pd.DataFrame(columns=['SNP position', 'gRNA', 'guide ID'])
    else:
        gRNA_df = pd.DataFrame(columns=['SNP position', 'gRNA', 'guide ID'])
    # we'll definitely want to include more of this information eventually. For now, let's just extract the sequence and the guideID
    snp_list=[]
    gRNA_list=[]
    id_list=[]
    for snp in postPAM_snps.pos:
        filt_gRNAs = gRNA_df[gRNA_df['SNP position']==snp]
        if len(filt_gRNAs)!=0:
            snp_list.append(snp)
            sequences=set(filt_gRNAs['gRNA'])
            gRNA_ids=set(filt_gRNAs['guide ID'])
            gRNA_list.append(sequences)
            id_list.append(gRNA_ids)
    final_df = pd.DataFrame({
        'pos':snp_list,
        'gRNA_sequence':gRNA_list,
        'gRNA_id':id_list
    })
    # merge back with data frame
    final_df = coerce_pos_numeric(final_df, 'pos')
    postPAM_snps = coerce_pos_numeric(postPAM_snps, 'pos')
    postPAM_snps = postPAM_snps.merge(final_df, on='pos',how='left')

    ## Finally lets merge pre- and post-PAM information

    # standardize data types for the merge
    postPAM_snps['CRISPR_Cas9_targetable']=True
    prePAM_snps = coerce_pos_numeric(prePAM_snps, 'pos')
    postPAM_snps = coerce_pos_numeric(postPAM_snps, 'pos')
    prePAM_snps['af'] = pd.to_numeric(prePAM_snps['af'], errors='coerce')
    postPAM_snps['af'] = pd.to_numeric(postPAM_snps['af'], errors='coerce')
    for col in ['ref', 'alt', 'edit_strat']:
        prePAM_snps[col] = prePAM_snps[col].astype('string')
        postPAM_snps[col] = postPAM_snps[col].astype('string')

    # merge
    all_snps = prePAM_snps.merge(postPAM_snps, on=['pos','af','ref','alt','edit_strat','heterozygote_freq',
                                                'het_num','homo_ref_freq','homo_ref_num','homo_alt_freq','homo_alt_num'],
                                how='outer')
    # fill the non targetable snps with 'False' in the last column
    all_snps['CRISPR_Cas9_targetable'] = (
        all_snps['CRISPR_Cas9_targetable']
        .where(all_snps['CRISPR_Cas9_targetable'].notna(), False)
        .astype(bool)
    )
    all_snps_reorder = all_snps[['pos','af','ref','alt','edit_strat','heterozygote_freq','homo_ref_freq','homo_alt_freq','het_num','homo_ref_num',
                                'homo_alt_num','CRISPR_Cas9_targetable','gRNA_sequence','gRNA_id']]

    # get current chromosome information
    exon_df=pd.read_csv(exon_file,index_col=0,dtype={"chromosome_name":"str"})
    cur_gene_chrom = exon_df[exon_df['hgnc_symbol']==gene]['chromosome_name'].values[0]
    # replace NAs in certain places with more descriptive logic
    all_snps_reorder['gRNA_sequence'] = all_snps_reorder['gRNA_sequence'].fillna("NA - No NGG PAMs nearby, not targetable by CRISPR/Cas9.")
    all_snps_reorder['gRNA_id'] = all_snps_reorder['gRNA_id'].fillna("NA")


    # Assign colors to snps
    maf_map={0.1:'106,109,133', 0.2:'0,0,0', 0.3:'0,0,0', 0.4:'0,26,255', 0.5:'0,26,255'} # grey, black, blue
    all_snps_reorder['MAF'] = np.minimum(all_snps_reorder['af'], 1 - all_snps_reorder['af'])
    colors=[]
    for idx,row in all_snps_reorder.iterrows():
        cur_maf = round(row.MAF,2)
        if cur_maf>=0.5:
            colors.append(maf_map[0.5])
        elif cur_maf>=0.4:
            colors.append(maf_map[0.4])
        elif cur_maf>=0.3:
            colors.append(maf_map[0.3])
        elif cur_maf>=0.2:
            colors.append(maf_map[0.2])
        elif cur_maf>=0.1:
            colors.append(maf_map[0.1])
        # else:
        #     print(cur_maf)
        #     print(round(cur_maf,2))
    all_snps_reorder['itemRgB']=colors

    # format the ref and alt alleles with the position
    all_snps_reorder['formatted_alleles'] = (
        'chr'
        + cur_gene_chrom
        + ":"
        + all_snps_reorder['pos'].astype(int).astype(str)
        + ' '
        + all_snps_reorder['ref']
        + '/'
        + all_snps_reorder['alt']
    )
    # format genotype frequencies
    all_snps_reorder['formatted_genos'] = all_snps_reorder.apply(
        lambda r:
            f"Heterozygous ({r.ref}{r.alt}) {int(r.heterozygote_freq * 100)}% ({r.het_num}/2548), "
            + f"Homozygous ref ({r.ref}{r.ref}) {int(r.homo_ref_freq * 100)}% ({r.homo_ref_num}/2548), "
            + f"Homozygous alt ({r.alt}{r.alt}) {int(r.homo_alt_freq * 100)}% ({r.homo_alt_num}/2548)",
        axis=1
    )

    
    ## Include population-specific frequencies
    # vcf formatting info
    num_samples=2548
    sample_list = list(range(1,num_samples+1)) # number of people in 1KG
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]

    def assert_unique_and_biallelic_vcf_values(vcf):
        # 1️⃣ remove fully duplicated rows
        vcf = vcf.drop_duplicates()
        
        # 2️⃣ keep only rows where ref and alt are single letters
        vcf = vcf[vcf['REF'].str.len() == 1]
        vcf = vcf[vcf['ALT'].str.len() == 1]
        
        # optional: reset index
        vcf = vcf.reset_index(drop=True)
        
        return vcf

    full_vcf = []
    for strat in strats:
        
        vcf_fp = os.path.join(results_dir, run_name, strat, "excavate/input_vcfs/" + gene + '_CommonVar_filtered.vcf.gz')
        if os.path.exists(vcf_fp) is False:
            continue
        cur_vcf = pd.read_table(vcf_fp, skiprows=23, compression='gzip')
        cur_vcf = assert_unique_and_biallelic_vcf_values(cur_vcf)
        full_vcf.append(cur_vcf)
    if len(full_vcf) > 0:
        full_vcf = pd.concat(full_vcf, ignore_index=True)
        full_vcf=assert_unique_and_biallelic_vcf_values(full_vcf)
    else:
        full_vcf = pd.DataFrame()
    # extract frequencies
    def parse_info(info_str):
        result = {}
        for item in info_str.split(";"):
            if "=" in item:
                k, v = item.split("=", 1)
                result[k] = v
            else:
                result[item] = True
        return result

    if full_vcf.empty:
        info_df = pd.DataFrame({'pos': all_snps_reorder['pos']})
        sample_names = []
    else:
        info_df = full_vcf["INFO"].apply(parse_info).apply(pd.Series)
        info_df["pos"] = full_vcf["POS"]
        # get sample names
        sample_names = list(full_vcf.columns[9:])
    # have a map of which sample belongs to which super-population
    sample_map_df = pd.read_table(sample_map_fp)
    sample_map = dict(zip(sample_map_df['Sample name'], sample_map_df['Superpopulation code']))
    # filter dictionary accordingly
    filt_map = {k: v for k, v in sample_map.items() if k in sample_names}
    # now calculate
    pop_snp_dict={}
    for pop in ['AFR', 'AMR', 'EAS', 'EUR', 'EUR,AFR', 'SAS']:
        if full_vcf.empty:
            pop_snp_dict[pop] = pd.DataFrame({
                'pos'+pop:[],
                'af'+pop:[],
                'pop_size'+pop:[],
                'het_freq'+pop:[],
                'het_num'+pop:[],
                'homo_ref_freq'+pop:[],
                'homo_ref_num'+pop:[],
                'homo_alt_freq'+pop:[],
                'homo_alt_num'+pop:[]
            })
            continue
        
        # get sample names for the current pop
        fixed_cols = list(full_vcf.columns[0:9])
        pop_filt_dict = {k: v for k, v in filt_map.items() if v in pop}
        cur_sample_cols=list(pop_filt_dict.keys())
        pop_size=len(cur_sample_cols)
        fixed_cols.extend(cur_sample_cols)
        # filter just to samples in the vcf and ones of that belong to the current pop
        
        pop_vcf = full_vcf[fixed_cols]
        # change first 9 columns case for function use
        pop_vcf.columns = ['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + list(pop_vcf.columns[9:])
        # now run existing functions to caclulate numbers
        
        # transpose and clean the vcf
        vcf_t, key_df = transpose_and_clean(pop_vcf)
        
        # calculate genotype frequencies
        het_dict=make_list_of_hets(vcf_t)
        ref_hom_dict=make_list_of_ref_homs(vcf_t)
        alt_hom_dict=make_list_of_alt_homs(vcf_t)
        
        # set up variables for storage
        snp_list=[]
        af_list=[]
        het_freq=[]
        het_num=[]
        homo_ref_freq=[]
        homo_ref_num=[]
        homo_alt_freq=[]
        homo_alt_num=[]
        # combine the data into general frequencies
        for snp,het_list in het_dict.items():
            # skips EUR,AFR merged pop
            if pop+'_AF' not in info_df.columns:
                continue
            snp_list.append(snp)
            af_list.append(info_df[info_df['pos']==snp][pop+'_AF'].values[0])
            het_num.append(len(het_list))
            het_freq.append(len(het_list)/pop_size if pop_size > 0 else np.nan)
            homo_ref_num.append(len(ref_hom_dict[snp]))
            homo_ref_freq.append(len(ref_hom_dict[snp])/pop_size if pop_size > 0 else np.nan)
            homo_alt_num.append(len(alt_hom_dict[snp]))
            homo_alt_freq.append(len(alt_hom_dict[snp])/pop_size if pop_size > 0 else np.nan)
        cur_pop_df = pd.DataFrame({
            'pos'+pop:snp_list,
            'af'+pop:af_list,
            'pop_size'+pop:pop_size,
            'het_freq'+pop:het_freq,
            'het_num'+pop:het_num,
            'homo_ref_freq'+pop:homo_ref_freq,
            'homo_ref_num'+pop:homo_ref_num,
            'homo_alt_freq'+pop:homo_alt_freq,
            'homo_alt_num'+pop:homo_alt_num
        })
        # store these variables for each population
        pop_snp_dict[pop]=cur_pop_df
        #all_snps_reorder = all_snps_reorder.merge(cur_pop_df, on=['pos'])
        
    # format the data frame for adding new columns
    all_snps_reorder['gRNA_id'] = all_snps_reorder['gRNA_id'].apply(str)
    all_snps_reorder['gRNA_sequence'] = all_snps_reorder['gRNA_sequence'].apply(str)
    all_snps_reorder.drop_duplicates(inplace=True, ignore_index=True)
    all_snps_reorder.sort_values(by='pos', ignore_index=True)
    # shorten syntax
    asr=all_snps_reorder.copy()
    for pop, df in pop_snp_dict.items():
        sorted_df=df.sort_values(by=['pos'+pop], ignore_index=True)
        asr = pd.concat([asr, sorted_df], axis=1)
    asr = coerce_pos_numeric(asr, 'pos')
    info_df = coerce_pos_numeric(info_df, 'pos')
    asr = asr.merge(info_df, on='pos', how='left')
    # format for final df
    for pop in ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']:

        het_freq = asr[f'het_freq{pop}']
        het_num  = asr[f'het_num{pop}']
        hr_freq  = asr[f'homo_ref_freq{pop}']
        hr_num   = asr[f'homo_ref_num{pop}']
        ha_freq  = asr[f'homo_alt_freq{pop}']
        ha_num   = asr[f'homo_alt_num{pop}']
        pop_size = asr[f'pop_size{pop}']
        afs = asr[f'{pop}_AF']

        asr[f'formatted_superpops_{pop}'] = (
            "Allele frequency = " + afs.astype(str) + ", " +
            "Heterozygous (" + asr['ref'] + asr['alt'] + ") " + (het_freq * 100).astype(int).astype(str) + "% (" + het_num.astype(str) + "/" + pop_size.astype(str) + "), "
            + "Homozygous ref (" + asr['ref'] + asr['ref'] + ") " + (hr_freq * 100).astype(int).astype(str) + "% (" + hr_num.astype(str) + "/" + pop_size.astype(str) + "), "
            + "Homozygous alt (" + asr['alt'] + asr['alt'] + ") " + (ha_freq * 100).astype(int).astype(str) + "% (" + ha_num.astype(str) + "/" + pop_size.astype(str) + ")"
        )

    # Finally, get the rsIDs
    rsIDs=[]
    if os.path.exists(rsID_fp):
        data = None
        max_attempts = 10
        retry_delay_seconds = 2
        for attempt in range(1, max_attempts + 1):
            try:
                # Open the file in read-binary mode ('rb')
                with gzip.open(rsID_fp, 'rb') as f:
                    # Load the data using pickle.load
                    data = pickle.load(f)
                break
            except (OSError, EOFError, pickle.UnpicklingError, gzip.BadGzipFile) as e:
                if attempt == max_attempts:
                    print(f"WARNING: Failed to read rsID map after {max_attempts} attempts for chromosome {cur_gene_chrom}: {e}")
                else:
                    time.sleep(retry_delay_seconds * attempt)
        if data is None:
            rsIDs = [np.nan] * len(asr)
        else:
            for pos in asr.pos:
                try:
                    cur_result = data.get(int(pos), np.nan)
                except (TypeError, ValueError):
                    cur_result = np.nan
                if isinstance(cur_result, str) and ',' in cur_result:
                    # take the first rsID if multiple
                    cur_result = cur_result.split(',')[0]
                rsIDs.append(cur_result)
    else:
        print(f"WARNING: rsID mapping file not found for chromosome {cur_gene_chrom}: {rsID_fp}")
        rsIDs = [np.nan] * len(asr)
    asr['rsID']=rsIDs
    # format the alleles
    asr['MAF'] = np.minimum(asr['af'], 1 - asr['af'])
    asr['formatted_name'] = asr['rsID'].fillna('') + ' ' + asr['ref'] + '/' + asr['alt'] + ', ' + 'MAF=' + asr['MAF'].round(2).astype(str)

    # finally, convert gene editing strategy names over to their more general terms:
    mapping = {
        "excision": "excision",
        "CRISPRoff": "epigenetic silencing",
        "donor_base_edits": "splice site disruption",
        "acceptor_base_edits":"splice site disruption",
        "indels": "exon disruption"
    }

    asr['edit_strat_set'] = asr['edit_strat'].apply(literal_eval)

    def recode_strategy(s):
        # s is a set of strings
        return {mapping.get(item, item) for item in s}  # default to original if not in mapping

    asr['edit_strat_general'] = asr['edit_strat_set'].apply(recode_strategy)
    asr['edit_strat_general_str'] = asr['edit_strat_general'].apply(lambda x: ", ".join(sorted(x)))

    # get 25 nucleotides +/- on either side of the variant
    fasta_filepath = ref_genome_fasta # "/wynton/group/capra/data/hg38_fasta/2022-03-14/hg38.fa.gz"
    fasta_open = pysam.Fastafile(fasta_filepath) # open the fasta file for an individual sample
    # get the position of the variants
    flanking_seqs=[]
    for pos in asr['pos']:
        seq = fasta_open.fetch(
            region=str(cur_gene_chrom),
            start=max(0, int(pos) - 26),
            end=int(pos) + 25,
        )
        if len(seq) == 51:
            formatted_seq = seq[:25] + "_" + seq[25] + "_" + seq[26:]
        else:
            formatted_seq = seq
        flanking_seqs.append(formatted_seq)
    asr['flanking_seqs'] = flanking_seqs

    ## Final formatting
    all_snps_reorder = asr.copy()
    bed_format = pd.DataFrame({
        'chrom':['chr' + cur_gene_chrom]*len(all_snps_reorder['pos']),
        'chromStart':all_snps_reorder['pos']-1,
        'chromEnd':all_snps_reorder['pos'],
        'name': all_snps_reorder['formatted_name'],
        'score':[0] *len(all_snps_reorder['pos']),
        'strand':['.'] * len(all_snps_reorder['pos']),
        'thickStart':all_snps_reorder['pos']-1,
        'thickEnd':all_snps_reorder['pos'],
        'itemRgb':all_snps_reorder['itemRgB'],
        'Description': ['A therapeutically editable common variant for the DnD gene ' + gene] * len(all_snps_reorder['pos']),
        'Ref/Alt':all_snps_reorder['ref']+'/'+all_snps_reorder['alt'],
        'AF':all_snps_reorder['af'].round(2),
        'MAF': all_snps_reorder['MAF'].round(2),
        'Global Genotype frequencies':all_snps_reorder['formatted_genos'],
        'AFR_genos':all_snps_reorder['formatted_superpops_AFR'],
        'EUR_genos':all_snps_reorder['formatted_superpops_EUR'],
        'AMR_genos':all_snps_reorder['formatted_superpops_AMR'],
        'EAS_genos':all_snps_reorder['formatted_superpops_EAS'],
        'SAS_genos':all_snps_reorder['formatted_superpops_SAS'],
        'Targetable by the following editing strategies': all_snps_reorder['edit_strat_general_str'],
        '5 closest excision partners, if any':[','.join(str(i) for i in closest_partners[x]) if x in closest_partners else 'NA' for x in all_snps_reorder['pos']],
        'CRISPR/SpCas9 Targetable?': ['Yes' if x==True else 'No' for x in all_snps_reorder['CRISPR_Cas9_targetable']],
        'Flanking sequences' : all_snps_reorder['flanking_seqs'],
        'Other Cas targetability tools':'CRISPOR - https://crispor.gi.ucsc.edu/, CRISPick - https://portals.broadinstitute.org/gppx/crispick/public'
    })
    bed_format['chromStart'] = bed_format['chromStart'].astype(int)
    bed_format['chromEnd'] = bed_format['chromEnd'].astype(int)
    bed_format['thickStart'] = bed_format['thickStart'].astype(int)
    bed_format['thickEnd'] = bed_format['thickEnd'].astype(int)
    bed_format['itemRgb'] = bed_format['itemRgb'].str.replace('"','')
    bed_format.sort_values(by=['chrom','chromStart'],ascending=True,inplace=True,ignore_index=True)

    # save
    # make path if it doesn't already exists
    bt_dir=args.bt_dir
    new_directory_name=os.path.join(bt_dir, 'per_gene_files', gene)
    try:
        os.makedirs(new_directory_name, exist_ok=True)
        print(f"Directory '{new_directory_name}' created successfully or already exists.")
    except OSError as e:
        print(f"Error creating directory: {e}")

    bed_format.to_csv(os.path.join(new_directory_name, gene + '_snp_track_ng.bed'),header=False,index=False, sep='\t')
    
    # create the .as file as well, for bigBed conversion:
    as_text = """table commonVar_geneEdits
    "BED file with CRISPR gene editing information"
    (
        string  chrom;         "Reference sequence chromosome or scaffold"
        uint    chromStart;    "Start position of feature on chromosome"
        uint    chromEnd;      "End position of feature on chromosome"
        string  name;          "Variant identifier"
        uint    score;         "Score"
        char[1] strand;        "+ or - for strand"
        uint    thickStart;    "Coding region start"
        uint    thickEnd;      "Coding region end"
        uint    itemRgb;       "Color corresponding to edit strategy"
        string  description;   "Description"
        string  RefAlt;        "Ref/Alt allele"
        string  af;            "Allele frequency (AF)"
        string  maf;           "Minor allele frequency (MAF)"
        lstring genos;         "Global genotype frequencies in 1000 Genomes"
        lstring afr;           "AFR genotype frequencies (1000 Genomes)"
        lstring eur;           "EUR genotype frequencies (1000 Genomes)"
        lstring amr;           "AMR genotype frequencies (1000 Genomes)"
        lstring eas;           "EAS genotype frequencies (1000 Genomes)"
        lstring sas;           "SAS genotype frequencies (1000 Genomes)"
        string  editStrats;    "Targetable by the following CRISPR editing strategies"
        string  exParts;       "Five closest excision partners, if any"
        string  PAMtargetable; "CRISPR/SpCas9 Targetable?"
        lstring flanking;      "+/- 25 bp flanking variant"
        string  castool;       "Cas/gRNA generation tools"
    )
    """

    # write one .as file for bigBed conversion
    if os.path.exists(os.path.join(bt_dir, 'metadata', 'bed_col_descriptors.as'))==False:
        with open(os.path.join(bt_dir, 'metadata', 'bed_col_descriptors.as'), "w") as f:
            f.write(as_text)


if __name__=="__main__":
    main()
