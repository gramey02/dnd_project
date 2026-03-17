# imports
import argparse
import os
import sys
import numpy as np
import pandas as pd
import pickle
from datetime import datetime
from itertools import permutations, combinations, product
import itertools
DATE = datetime.now().strftime("%Y-%m-%d")
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from collections import defaultdict

# create haplotype df where each person's two haplotypes are two rows
def create_haplotype_df(vcf, string_list):
    vcf=vcf.copy()
    if 'edit_strategy' in vcf.columns:
        vcf.drop(labels=['edit_strategy'], axis=1, inplace=True)
    
    # more efficient haplotype assignment
    people = ["sample" + s for s in string_list]
    
    sample_haps_suffix = (
        [f"{p}.1" for p in people] +
        [f"{p}.2" for p in people]
    )
    sample_haps_suffix.sort()
    
    geno = vcf[people]  # shape: (n_sites, n_people)
    
    hap1 = geno.apply(lambda col: col.str.split("|", expand=True)[0])
    hap2 = geno.apply(lambda col: col.str.split("|", expand=True)[1])
    
    hap1 = hap1.astype(float)
    hap2 = hap2.astype(float)
    
    # homozygous → NaN
    mask_hom = hap1 == hap2
    hap1[mask_hom] = np.nan
    hap2[mask_hom] = np.nan
    
    pos_cols = vcf["pos"].astype(str).values
    
    df = pd.concat(
        [
            hap1.T.rename(index=lambda x: f"{x}.1"),
            hap2.T.rename(index=lambda x: f"{x}.2"),
        ],
        axis=0
    )
    
    df.columns = pos_cols
    # remove duplicate columns
    df = df.loc[:, ~df.columns.duplicated()].copy()

    # extract numeric person id and haplotype as integers
    tmp = df.index.to_series().str.extract(r'sample(\d+)\.(\d)')
    person_id = tmp[0].astype(int)
    hap_id = tmp[1].astype(int)
    
    # sort by person id, then haplotype
    df = df.iloc[np.lexsort((hap_id, person_id))]

    # run check to make sure heterozygotes were properly stored in df
    for col in df.columns:
        # get het count for the snp pos in the vcf
        first_snp = vcf.loc[(vcf["pos"] == int(col))].iloc[0]
    
        vcf_het_count = first_snp.loc[
            first_snp.index.str.startswith("sample")
        ].isin(["0|1", "1|0"]).sum()
    
        # get het count for the snp pos in the df
        df_het_count = (df.loc[:, col].notna().sum())/2
        
        if df_het_count!=vcf_het_count:
            raise ValueError(
            f"Count mismatch (SNP: {col}): VCF_het_count={vcf_count}, DF_het_count={df_count}"
        )

    return df

def base_sample(idx):
    return idx.rsplit(".", 1)[0]

def iterative_pick_and_prune_with_tracking(df, gene, output_dir, strat):
    # check if a checkpoint already exists
    ckpt_dir=os.path.join(output_dir, 'summary_files', 'cross_strat_gRNAs', 'checkpoints', gene + '_'+strat+'_checkpoint.pkl')
    if os.path.exists(ckpt_dir):
        picks, rows_removed_per_iter, samples_removed_per_iter, df, remaining_by_sample, iteration = load_checkpoint_non_excision(ckpt_dir)
    # otherwise start from the inputted df
    else:
        df = df.copy()
        picks = []
        rows_removed_per_iter = []
        samples_removed_per_iter = []
        # Track remaining rows per base sample
        remaining_by_sample = defaultdict(set)
        for idx in df.index:
            remaining_by_sample[base_sample(idx)].add(idx)
        iteration=0
    
    while True:
        best = None
        best_key = None

        for col_idx, col in enumerate(df.columns):
            counts = df[col].value_counts(dropna=True)

            for val in (0.0, 1.0):
                if val not in counts:
                    continue

                count = counts[val]
                key = (count, val == 0.0, -col_idx)

                if best_key is None or key > best_key:
                    best_key = key
                    best = (count, col, val)

        # Stop when only NaNs remain
        if best is None:
            break

        count, col, val = best
        picks.append((col, val, count))

        # Identify rows to remove this iteration
        to_remove = df.index[df[col] == val].tolist()
        rows_removed_per_iter.append(len(to_remove))

        fully_removed_samples = set()

        for idx in to_remove:
            base = base_sample(idx)
            remaining_by_sample[base].discard(idx)

            # If both .1 and .2 are now gone → full sample removed
            if len(remaining_by_sample[base]) == 0:
                fully_removed_samples.add(base)

        samples_removed_per_iter.append(sorted(fully_removed_samples))

        # Remove rows
        df = df.drop(index=to_remove)

        # save checkpoint
        if iteration%5==0:
            create_checkpoint_non_excision(picks, rows_removed_per_iter, samples_removed_per_iter, df, remaining_by_sample, ckpt_dir, iteration)
        iteration+=1

    return {
        "final_df": df,
        "picks": pd.DataFrame(picks, columns=['snp_pos','ref_alt','num_haplotypes_targeted']),
        "rows_removed_per_iteration": rows_removed_per_iter,
        "samples_removed_per_iteration": samples_removed_per_iter,
    }

def reformat_non_excision_summary_df(result):
    # convert the format of this summary df
    summary_df = result['picks']
    snp_alleles=[]
    for idx,row in summary_df.iterrows():
        snp_alleles.append(str(row.snp_pos) + "_" + str(int(row.ref_alt)))
    summary_df = summary_df.copy()
    summary_df['snp1_allele'] = snp_alleles
    summary_df['snp2_allele'] = None
    
    # get the cumulative number of hets targeted at both haplotypes upon each iteration
    cumulative_hets=0
    cumulative_hets_list=[]
    for item in result['samples_removed_per_iteration']:
        cumulative_hets=cumulative_hets+len(item)
        cumulative_hets_list.append(cumulative_hets)
    
    # get the cumulative number of haplotypes targeted upon each iteration
    cumulative_haps=0
    cumulative_haps_list=[]
    for item in result['picks']['num_haplotypes_targeted']:
        cumulative_haps=cumulative_haps+item
        cumulative_haps_list.append(cumulative_haps)
    
    summary_df['num_guides']=list(range(1,len(summary_df.snp_pos)+1))
    summary_df['num_hets_selected']=[len(x) for x in result['samples_removed_per_iteration']]
    summary_df['cumulative_hets_fully_targeted'] = cumulative_hets_list
    summary_df['cumulative_haplotypes_targeted'] = cumulative_haps_list
    
    summary_df = summary_df[['snp1_allele','snp2_allele','num_haplotypes_targeted','num_hets_selected','cumulative_haplotypes_targeted','cumulative_hets_fully_targeted', 'num_guides']]
    summary_df.rename(columns={'num_haplotypes_targeted':'haplotypes_added','num_hets_selected':'people_added', 'cumulative_hets_fully_targeted':'cumulative_people_targeted'},inplace=True)
    
    return summary_df

def create_checkpoint_non_excision(picks, rows_removed_per_iter, samples_removed_per_iter, df, remaining_by_sample, ckpt_dir, iteration):
    ckpt=[picks, rows_removed_per_iter, samples_removed_per_iter, df, remaining_by_sample, iteration]
    with open(ckpt_dir, 'wb') as fp:
        pickle.dump(ckpt, fp)
    return None

def load_checkpoint_non_excision(ckpt_fp):
    with open(ckpt_fp, 'rb') as fp:
        ckpt = pickle.load(fp)
    return ckpt[0], ckpt[1], ckpt[2], ckpt[3], ckpt[4], ckpt[5]

def assert_unique_and_biallelic_vcf_values(vcf):
    # 1️⃣ remove fully duplicated rows
    vcf = vcf.drop_duplicates()
    
    # 2️⃣ keep only rows where ref and alt are single letters
    vcf = vcf[vcf['ref'].str.len() == 1]
    vcf = vcf[vcf['alt'].str.len() == 1]
    
    # optional: reset index
    vcf = vcf.reset_index(drop=True)

    return vcf

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type = str, required = True, help = 'Output directory for het information.')
    parser.add_argument('--gene', type = str, required = True, help = 'Current gene to prioritize guides for.')
    parser.add_argument('--num_samples', type=int, required=True, help = 'Number of samples in vcf files')
    parser.add_argument('--all_strats_together', type=str, required=True, help='Run all non-excision strats together or separately.')
    args = parser.parse_args()
    return args

def main():
    args=parse_args()
    strats=['indels','CRISPRoff','donor_base_edits','acceptor_base_edits']
    num_samples=args.num_samples
    output_dir = args.output_dir
    gene=args.gene
    all_strats_together = args.all_strats_together

    # set up variables for loading vcf files
    sample_list = list(range(1,num_samples+1)) # number of people in 1KG
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]

    # load the gene's targetable snps for each edit strat
    vcf=None

    # run all non-excision strats together or separately
    if (all_strats_together==1) or (all_strats_together=='True') or (all_strats_together==True):
        for strat in strats:
            filtered_vcf_dir = os.path.join(output_dir, strat, 'excavate/Guide_filtered_vcfs')
            fp_to_check = os.path.join(output_dir, strat, 'excavate/Guide_filtered_vcfs', gene + '_guide_filtered.vcf')
            if (gene + '_guide_filtered.vcf') in os.listdir(filtered_vcf_dir):
                # load gene's snps
                cur_vcf = pd.read_table(fp_to_check, comment='#', header=None)
                cur_vcf.columns=cols
                cur_vcf['edit_strategy']=strat
                vcf=pd.concat([vcf, cur_vcf])

        if vcf is None or vcf.empty:
            extra="No vcf loaded for all_non_excision_strats."
            log_dir = os.path.join(output_dir, "summary_files/cross_strat_gRNAs/logs")
            log_fp = os.path.join(log_dir, f"{gene}.log")
            if not os.path.exists(log_fp):
                with open(log_fp, "w") as f:
                    f.write("timestamp\tgene\tremaining_haps\textra\n")
            os.makedirs(log_dir, exist_ok=True)
            log_fp = os.path.join(log_dir, f"{gene}.log")
            ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            line = (
                f"{ts}\t"
                f"{gene}\t"
            )
            if extra:
                line += f"\t{extra}"
            with open(log_fp, "a") as f:
                f.write(line + "\n")

        else:
            vcf = assert_unique_and_biallelic_vcf_values(vcf)
            # truncate vcf just to make sure a smaller version of this code runs okay
            # vcf=vcf.copy()
            # vcf=vcf[0:15]
            #print('finished truncating vcf file')

            # create haplotype-based df
            df = create_haplotype_df(vcf, string_list)
            print('Finished creating data frame for ' + gene)

            result = iterative_pick_and_prune_with_tracking(df, gene, output_dir, 'all_strats')
            summary_df = reformat_non_excision_summary_df(result)

            summary_df.to_csv(os.path.join(output_dir, 'summary_files', 'cross_strat_gRNAs', 'results',gene+'_non_excision_all_strats_gRNAs.csv'))
            # save the resulting summary df below
            if summary_df is not None:
                # note if the job completed or failed
                success_fp = os.path.join(
                    output_dir,
                    "summary_files",
                    "cross_strat_gRNAs",
                    "logs",
                    f"{gene}.DONE_all_non_excision"
                )
                with open(success_fp, "w") as f:
                    f.write("OK\n")
    else:
        for strat in strats:
            vcf=None
            filtered_vcf_dir = os.path.join(output_dir, strat, 'excavate/Guide_filtered_vcfs')
            fp_to_check = os.path.join(output_dir, strat, 'excavate/Guide_filtered_vcfs', gene + '_guide_filtered.vcf')
            if (gene + '_guide_filtered.vcf') in os.listdir(filtered_vcf_dir):
                # load gene's snps
                cur_vcf = pd.read_table(fp_to_check, comment='#', header=None)
                cur_vcf.columns=cols
                cur_vcf['edit_strategy']=strat
                vcf=cur_vcf.copy()

            # if no vcf is loaded, note it in a log
            if vcf is None or vcf.empty:
                continue

            vcf = assert_unique_and_biallelic_vcf_values(vcf)
            # truncate vcf just to make sure a smaller version of this code runs okay
            # vcf=vcf.copy()
            # vcf=vcf[0:15]
            #print('finished truncating vcf file')

            # create haplotype-based df
            df = create_haplotype_df(vcf, string_list)
            print('Finished creating data frame for ' + gene)

            result = iterative_pick_and_prune_with_tracking(df, gene, output_dir, strat)
            summary_df = reformat_non_excision_summary_df(result)

            summary_df.to_csv(os.path.join(output_dir, 'summary_files', 'cross_strat_gRNAs', 'results',gene+'_' + strat + '_gRNAs.csv'))    
            # note if it finished running
            if summary_df is not None:
                # note if the job completed or failed
                success_fp = os.path.join(
                    output_dir,
                    "summary_files",
                    "cross_strat_gRNAs",
                    "logs",
                    f"{gene}.DONE_{strat}"
                )
                with open(success_fp, "w") as f:
                    f.write("OK\n")

if __name__=='__main__':
    main()