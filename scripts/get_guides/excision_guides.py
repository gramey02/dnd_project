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
import math

# lots of helper functions

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
            f"Count mismatch (SNP: {col}): VCF_het_count={vcf_het_count}, DF_het_count={df_het_count}"
        )

    return df

def df_to_numpy_with_index(df):
    """
    Convert a DataFrame to a dict of NumPy arrays per column and a NumPy array of row indices.

    Returns:
        col_arrays: dict mapping column name -> np.array of values
        h_index: np.array of original row indices
    """
    col_arrays = {col: df[col].to_numpy() for col in df.columns}
    h_index = df.index.to_numpy()
    return col_arrays, h_index


# double snp version
def iterative_pick_and_prune_pairs_with_tracking(df, valid_pairs):
    df = df.copy()

    # Normalize valid pairs
    valid_pair_set = {
        tuple(sorted(map(str, pair)))
        for pair in valid_pairs
    }

    # Keep only pairs that are actually present in df
    columns_set = set(df.columns)

    pairs_to_check = [
        (a, b)
        for (a, b) in valid_pair_set
        if a in columns_set and b in columns_set
    ]

    # Cache column order lookup (avoid repeated get_loc)
    col_positions = {col: i for i, col in enumerate(df.columns)}

    picks = []
    rows_removed_per_iter = []
    samples_removed_per_iter = []

    def base_sample(idx):
        return idx.rsplit(".", 1)[0]

    # Precompute sample tracking
    remaining_by_sample = defaultdict(set)
    for idx in df.index:
        remaining_by_sample[base_sample(idx)].add(idx)

    while True:
        best = None
        best_key = None

        columns = list(df.columns)

        # Iterate over column pairs
        for col_i, col_j in pairs_to_check:

            # Enforce valid SNP pairing
            if tuple(sorted((str(col_i), str(col_j)))) not in valid_pair_set:
                continue

            # Pull NumPy arrays once (FAST)
            col1 = df[col_i].values
            col2 = df[col_j].values

            # Mask for non-NaN rows (replaces dropna)
            valid_mask = (~np.isnan(col1)) & (~np.isnan(col2))
            if not valid_mask.any():
                continue

            c1 = col1[valid_mask]
            c2 = col2[valid_mask]

            # Count haplotypes directly (avoid inner product loop)
            counts = [
                (np.sum((c1 == 0.0) & (c2 == 0.0)), (0.0, 0.0)),
                (np.sum((c1 == 0.0) & (c2 == 1.0)), (0.0, 1.0)),
                (np.sum((c1 == 1.0) & (c2 == 0.0)), (1.0, 0.0)),
                (np.sum((c1 == 1.0) & (c2 == 1.0)), (1.0, 1.0)),
            ]

            for count, (v_i, v_j) in counts:
                if count == 0:
                    continue

                key = (
                    count,                                   # maximize removed rows
                    (v_i == 0.0) + (v_j == 0.0),             # prefer more 0s
                    -col_positions[col_i],                   # deterministic tie-break
                    -col_positions[col_j],
                )

                if best_key is None or key > best_key:
                    best_key = key
                    best = (count, col_i, col_j, (v_i, v_j))

        # Stop when no valid pattern found
        if best is None:
            break

        count, col_i, col_j, (v_i, v_j) = best
        picks.append((col_i, col_j, (v_i, v_j), count))

        # Use NumPy mask again (no new Series objects)
        col1 = df[col_i].values
        col2 = df[col_j].values
        remove_mask = (col1 == v_i) & (col2 == v_j)

        to_remove = df.index[remove_mask]

        rows_removed_per_iter.append(len(to_remove))

        fully_removed_samples = set()
        for idx in to_remove:
            base = base_sample(idx)
            remaining_by_sample[base].discard(idx)
            if len(remaining_by_sample[base]) == 0:
                fully_removed_samples.add(base)

        samples_removed_per_iter.append(sorted(fully_removed_samples))

        # Drop rows once
        df = df.drop(index=to_remove)

        # You break after one iteration in original code
        break

    return {
        "final_df": df,
        "picks": pd.DataFrame(
            picks,
            columns=["snp1", "snp2", "haplotype", "num_haplotypes_targeted"]
        ),
        "rows_removed_per_iteration": rows_removed_per_iter,
        "samples_removed_per_iteration": samples_removed_per_iter,
    }

def pick_best_per_allele_yield(sg_best_allele, sg_best_haplo_num, sg_best_haplos, dg_best_allele_pair, dg_best_haplo_num, dg_best_haplos):
    # Handle missing cases
    if sg_best_allele is None and dg_best_allele_pair is None:
        return set(), (), 0

    if dg_best_allele_pair is None:
        return sg_best_haplos, sg_best_allele, sg_best_haplo_num

    if sg_best_allele is None:
        return dg_best_haplos, dg_best_allele_pair, dg_best_haplo_num

    # Compare per-allele yield
    if dg_best_haplo_num / 2 > sg_best_haplo_num:          # note that this forces ties to favor the single guide case
        return dg_best_haplos, dg_best_allele_pair, dg_best_haplo_num
    else:
        return sg_best_haplos, sg_best_allele, sg_best_haplo_num
    
def get_untargeted_alleles(col_arrays, h_index, targeted_alleles):
    """
    Args:
        col_arrays: dict of {snp_col_name: np.array of haplotype values (float)}
        h_index: np.array of row indices (haplotype names)
        targeted_alleles: list of strings in the form "snp_allele"
    
    Returns:
        List of untargeted allele strings, e.g., ["5_0", "5_1", "7_1"]
    """
    untargeted_alleles=[]
    for snp, col in col_arrays.items():
        # get unique non-NaN alleles in the column
        unique_vals = np.unique(col[~np.isnan(col)])
        alleles = [f"{snp}_{int(val)}" for val in unique_vals]
        untargeted_alleles.extend(alleles)
    untargeted_alleles_final_filt = [x for x in untargeted_alleles if x not in targeted_alleles]
    return untargeted_alleles_final_filt

def numeric_key(s):
    if not isinstance(s, str) or "_" not in s:
        raise ValueError(f"Invalid allele format: {s}")
    return int(s.split("_", 1)[0])
def allele_key(s):
    return int(s.split("_", 1)[1])

def filter_to_valid_pairs(allele_pairs, valid_pairs):
    filtered = []

    for pair in allele_pairs:
        a = numeric_key(pair[0])
        b = numeric_key(pair[1])

        if a <= b:
            snp_pair = (a, b)
        else:
            snp_pair = (b, a)

        if snp_pair in valid_pairs:
            filtered.append(pair)  # keep ORIGINAL allele strings

    return filtered
    
def get_targeted_haplos(targeted_alleles, a, col_arrays, h_index, valid_pairs):
    # if allele is a string, it means there is only one allele as input
    if isinstance(a, str):
        # make all pairs of one previously targeted allele and new allele
        allele_pairs = [tuple(sorted((x, a), key=numeric_key)) for x in targeted_alleles]

    # if allele is tuple, it means there are two new alleles as input
    if isinstance(a,tuple):
        # make all pairs of one previously targeted allele and each of the two new alleles
        allele_pairs = [tuple(sorted((x, y), key=numeric_key)) for x in targeted_alleles for y in a]
        # also add the inputted tuple itself as an additional pair
        allele_pairs = allele_pairs + [a]
        
    # deduplicate any pairs
    allele_pairs = list(set(allele_pairs))
    # filter resulting pairs to valid pairs
    filtered_allele_pairs=filter_to_valid_pairs(allele_pairs, valid_pairs)

    # get targetable haplotypes for each allele pair
    new_targetable_haplos=set()
    for pair in filtered_allele_pairs:
        # get the snps and alleles in the pair
        snp1, snp2 = str(numeric_key(pair[0])), str(numeric_key(pair[1]))
        a1, a2 = allele_key(pair[0]), allele_key(pair[1])
        
        # filter the columns and rows of the haplotype df to get the people targetable by the current pair
        col1 = col_arrays[snp1]
        col2 = col_arrays[snp2]

        # mask rows that match allele pair and are non-NaN
        mask = (~np.isnan(col1)) & (~np.isnan(col2)) & (col1 == a1) & (col2 == a2)
        new_targetable_haplos.update(h_index[mask])
        
    # return the union of targetable_hets
    return new_targetable_haplos

def single_guide_addition(untargeted_alleles, col_arrays, h_index, targeted_alleles, valid_pairs):

    # account for nans
    if not untargeted_alleles:
        return None, 0, set()

    
    new_targeted_haplos_len_per_added_allele=[]
    new_targeted_haplos_per_added_allele=[]
    selected_alleles=[]

    # loop through untargeted alleles, adding one at a time to existing pairs, and see how many additional haplotypes are captured by adding it
    for a in untargeted_alleles:
        new_targetable_haplos = get_targeted_haplos(targeted_alleles, a, col_arrays, h_index, valid_pairs)
        selected_alleles.append(a)
        new_targeted_haplos_len_per_added_allele.append(len(new_targetable_haplos))
        new_targeted_haplos_per_added_allele.append(new_targetable_haplos)

    if not new_targeted_haplos_len_per_added_allele:
        return None, 0, set()
        
    # select the added allele that targets the most additional haplos
    max_idx, max_val = max(
        enumerate(new_targeted_haplos_len_per_added_allele),
        key=lambda x: x[1]
    )
    
    best_allele = selected_alleles[max_idx]
    best_haplo_num = max_val
    best_haplos = new_targeted_haplos_per_added_allele[max_idx]

    
    return best_allele, best_haplo_num, best_haplos

def dual_guide_addition(untargeted_alleles, col_arrays, h_index, targeted_alleles, valid_pairs):
    new_targeted_haplos_len_per_added_pair=[]
    new_targeted_haplos_per_added_pair=[]
    selected_allele_pairs=[]

    # guard against case where less than 2 untargeted alleles remain
    if len(untargeted_alleles) < 2:
        return None, 0, set()

    # loop through untargeted alleles, adding one pair at a time to existing pairs, and see how many additional haplotypes are captured by adding the two
    pairs = [tuple(sorted(pair, key=numeric_key)) for pair in combinations(untargeted_alleles, 2)]
    for a_pair in pairs:
        if (numeric_key(a_pair[0]), numeric_key(a_pair[1])) not in valid_pairs:       # checking if the pair is even valid
            continue
        new_targetable_haplos = get_targeted_haplos(targeted_alleles, a_pair, col_arrays, h_index, valid_pairs)
        selected_allele_pairs.append(a_pair)
        new_targeted_haplos_len_per_added_pair.append(len(new_targetable_haplos))
        new_targeted_haplos_per_added_pair.append(new_targetable_haplos)

    if not new_targeted_haplos_len_per_added_pair:
        return None, 0, set()

    # select the added pair that targets the most additional haplos
    max_idx, max_val = max(
        enumerate(new_targeted_haplos_len_per_added_pair),
        key=lambda x: x[1]
    )

    best_allele_pair = selected_allele_pairs[max_idx]
    best_haplo_num = max_val
    best_haplos = new_targeted_haplos_per_added_pair[max_idx]

    return best_allele_pair, best_haplo_num, best_haplos

# gets sample name without its haplotype appended
def base_sample(idx):
    return idx.rsplit(".", 1)[0]

def create_summary_df(allele_additions_per_iter, haplotype_additions_per_iter, sample_additions_per_iter):
    sum_df = pd.DataFrame({
        "snp1_allele": pd.Series([None] * len(allele_additions_per_iter), dtype="object"),
        "snp2_allele": pd.Series([None] * len(allele_additions_per_iter), dtype="object"),
        "haplotypes_added": pd.Series([np.nan] * len(allele_additions_per_iter), dtype="Int64"),
        "people_added": pd.Series([np.nan] * len(allele_additions_per_iter), dtype="Int64"),
    })
    haplotypes_added = [len(y) for y in haplotype_additions_per_iter]
    people_added = [len(x) for x in sample_additions_per_iter]
    
    num_guides_list=[]
    guide_counter=0
    for idx,row in sum_df.iterrows():
        # fill in snp information
        if isinstance(allele_additions_per_iter[idx], tuple):
            sum_df.loc[idx,'snp1_allele'] = allele_additions_per_iter[idx][0]
            sum_df.loc[idx,'snp2_allele'] = allele_additions_per_iter[idx][1]
            guide_counter+=2
            num_guides_list.append(guide_counter)
        else:
            sum_df.loc[idx,'snp1_allele'] = allele_additions_per_iter[idx]
            guide_counter+=1
            num_guides_list.append(guide_counter)
    
        # fill in haplotype and sample information
        sum_df.loc[idx,'haplotypes_added'] = haplotypes_added[idx]
        sum_df.loc[idx,'people_added'] = people_added[idx]
    
    # combine some information to get cumulative numbers
    
    # get the cumulative number of haplotypes targeted upon each iteration
    cumulative_haps=0
    cumulative_haps_list=[]
    for item in haplotypes_added:
        cumulative_haps=cumulative_haps+item
        cumulative_haps_list.append(cumulative_haps)
    
    # get the cumulative number of people added upon each iteration
    cumulative_hets=0
    cumulative_hets_list=[]
    for item in people_added:
        cumulative_hets=cumulative_hets+item
        cumulative_hets_list.append(cumulative_hets)
    
    sum_df['cumulative_haplotypes_targeted'] = cumulative_haps_list
    sum_df['cumulative_people_targeted'] = cumulative_hets_list
    sum_df['num_guides'] = num_guides_list
    return sum_df

def assert_unique_and_biallelic_vcf_values(vcf):
    # 1️⃣ remove fully duplicated rows
    vcf = vcf.drop_duplicates()
    
    # 2️⃣ keep only rows where ref and alt are single letters
    vcf = vcf[vcf['ref'].str.len() == 1]
    vcf = vcf[vcf['alt'].str.len() == 1]
    
    # optional: reset index
    vcf = vcf.reset_index(drop=True)

    return vcf

def load_checkpoint(checkpoint_dir):
    with open(checkpoint_dir, 'rb') as fp:
        ckpt_data = pickle.load(fp)
    return ckpt_data[0], ckpt_data[1], ckpt_data[2], ckpt_data[3], ckpt_data[4], ckpt_data[5], ckpt_data[6], ckpt_data[7]

def create_checkpoint(allele_additions_per_iter, haplotype_additions_per_iter, sample_additions_per_iter, h_index, col_arrays, remaining_by_sample, targeted_alleles, iteration, gene, output_dir):
    ckpt_list=[allele_additions_per_iter, haplotype_additions_per_iter, sample_additions_per_iter, h_index, col_arrays, remaining_by_sample, targeted_alleles, iteration]
    checkpoint_dir=output_dir + '/summary_files/cross_strat_gRNAs/excision_guides/checkpoints/' + gene + '_ckpt.pkl'
    with open(checkpoint_dir, 'wb') as fp:
        pickle.dump(ckpt_list, fp)
    return None

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type = str, required = True, help = 'Output directory for het information.')
    parser.add_argument('--gene', type = str, required = True, help = 'Current gene to prioritize guides for.')
    parser.add_argument('--num_samples', type=int, required=True, help = 'Number of samples in vcf files')
    parser.add_argument('--valid_pairs_fp', type=str, required=True, help='Run all non-excision strats together or separately.')
    parser.add_argument('--vcf_dir', type=str, required=True, help='VCF directory for guide-filtered vcfs.')
    parser.add_argument('--max_iter', type=int, required=True, help='max number of iterations on which to pick gRNAs.')
    args = parser.parse_args()
    return args

def main():
    args=parse_args()
    num_samples=args.num_samples
    output_dir = args.output_dir
    gene=args.gene
    valid_pairs_fp = args.valid_pairs_fp
    valid_pairs_fp=valid_pairs_fp+'/'+gene+'_valid_snp_pairs.pkl'
    vcf_dir=args.vcf_dir
    max_iter=args.max_iter

    # gene='ABCC8' # 'NEFL'
    # num_samples=2548
    # output_dir = "/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/pipeline_results/RUN_12_19_25"
    # valid_pairs_fp = "/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/pipeline_results/RUN_12_19_25/excision/CommonVars/valid_snp_pairs"
    # valid_pairs_fp=valid_pairs_fp+'/'+gene+'_valid_snp_pairs.pkl'
    # vcf_dir="/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/pipeline_results/RUN_12_19_25/excision/excavate/Guide_filtered_vcfs"
    # max_iter=50

    # set up variables for loading vcf files
    sample_list = list(range(1,num_samples+1)) # number of people in 1KG
    string_list = list(map(str, sample_list))
    cols=['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format'] + ["sample" + s for s in string_list]

    # load the gene's targetable snps for each edit strat
    vcf = pd.read_table(vcf_dir + "/" + gene + '_guide_filtered.vcf', comment='#', header=None) # pd.read_table(vcf_dir + "/" + gene + '_guide_filtered.vcf', comment='#', header=None)
    vcf.columns=cols
    vcf=assert_unique_and_biallelic_vcf_values(vcf)
    print('loaded vcf')

    # load valid snp pairs for the current gene
    with open(valid_pairs_fp, 'rb') as fp:
        valid_pairs = pickle.load(fp)
    valid_pairs=set(valid_pairs) # this allows for faster lookup in the subsequent functions

    # load a checkpoint if it already exists
    ckpt_dir=output_dir + '/summary_files/cross_strat_gRNAs/excision_guides/checkpoints/' + gene + '_ckpt.pkl'
    if os.path.exists(ckpt_dir):
        allele_additions_per_iter, haplotype_additions_per_iter, sample_additions_per_iter, h_index, col_arrays, remaining_by_sample, targeted_alleles, iteration = load_checkpoint(ckpt_dir)
    else:
        df=create_haplotype_df(vcf, string_list)
        # run algorithm to pick the first pair of alleles
        haplotype_names=set(list(df.index))
        result = iterative_pick_and_prune_pairs_with_tracking(df, valid_pairs)

        # store these picked alleles
        snp1 = result['picks']['snp1'].values[0]
        snp2 = result['picks']['snp2'].values[0]
        a1,a2 = result['picks']['haplotype'].values[0]
        remaining_haps = result['final_df']
        col_arrays, h_index = df_to_numpy_with_index(remaining_haps)

        # get the first two targeted alleles
        targeted_alleles = [f"{snp1}_{int(a1)}", f"{snp2}_{int(a2)}"]

        # tracks the allele(s) added at each iteration
        allele_additions_per_iter=[]
        allele_additions_per_iter.append((f"{snp1}_{int(a1)}", f"{snp2}_{int(a2)}"))

        # tracks the haplotypes added at each iteration
        haplotype_additions_per_iter=[]
        haplotype_additions_per_iter.append(haplotype_names - set(h_index)) # add those added from first pair here

        # tracks the samples (aka people) captured at each iteration
        remaining_by_sample = defaultdict(set)
        for idx in haplotype_names:
            remaining_by_sample[base_sample(idx)].add(idx)
        sample_additions_per_iter = []
        sample_additions_per_iter.append(result['samples_removed_per_iteration'][0]) # add those captured from the first pair here (should be 0, since only one haplotype should be captured so far)
        # need to remove the already picked haplos from remaining_by_sample
        for idx in (haplotype_names - set(h_index)):
            base = base_sample(idx)
            if idx in remaining_by_sample[base]:
                remaining_by_sample[base].remove(idx)
        iteration=0

    #print('first pair selected')

    while len(h_index)>0 and (iteration<max_iter):
        # get number of remaining haplotypes
        prev_remaining = len(h_index)

        # get untargeted alleles remaining
        untargeted_alleles = get_untargeted_alleles(col_arrays, h_index, targeted_alleles) # this should format and remove the already two targeted alleles from the set
        
        # calculate the benefit from adding a single guide to the pair
        sg_best_allele, sg_best_haplo_num, sg_best_haplos = single_guide_addition(untargeted_alleles, col_arrays, h_index, targeted_alleles, valid_pairs)
        
        # calculate the benefit from adding a completely new pair of guides
        dg_best_allele_pair, dg_best_haplo_num, dg_best_haplos = dual_guide_addition(untargeted_alleles, col_arrays, h_index, targeted_alleles, valid_pairs)
        
        # pick the best per-allele yield of haplotypes
        best_haplos=None
        best_allele_s = None
        best_haplo_num = None
        best_haplos, best_allele_s, best_haplo_num = pick_best_per_allele_yield(sg_best_allele, sg_best_haplo_num, sg_best_haplos, dg_best_allele_pair, dg_best_haplo_num, dg_best_haplos)
        
        # break poinnt for if no allele was selected
        if best_allele_s is None:
            break
        # another catch-all break point: if no new haplotypes captured, then break
        if not best_haplos:
            break
            
        # update selections and haplotype data frame accordingly
        allele_additions_per_iter.append(best_allele_s)
        best_allele_list = list(best_allele_s) if isinstance(best_allele_s, tuple) else [best_allele_s]
        targeted_alleles.extend(best_allele_list)
        haplotype_additions_per_iter.append(best_haplos)
        # Track newly fully captured samples this iteration
        newly_fully_removed = set()
        for idx in best_haplos:
            base = base_sample(idx)
            if idx in remaining_by_sample[base]:
                remaining_by_sample[base].remove(idx)
                # Sample is fully targeted only when *all* its haplotypes are gone
                if len(remaining_by_sample[base]) == 0:
                    newly_fully_removed.add(base)
        # Record per-iteration sample additions
        sample_additions_per_iter.append(sorted(newly_fully_removed))
        
        # --- remove captured haplos from NumPy arrays ---
        keep_mask = ~np.isin(h_index, list(best_haplos))
        h_index = h_index[keep_mask]
        for col in col_arrays:
            col_arrays[col] = col_arrays[col][keep_mask]

        # store a checkpoint
        if iteration%5==0:
            create_checkpoint(allele_additions_per_iter, haplotype_additions_per_iter, sample_additions_per_iter, h_index, col_arrays, remaining_by_sample, targeted_alleles, iteration, gene, output_dir)
        # update iteration
        iteration+=1

        # another catch-all break point: if no haplotypes were removed from the data frame on this iteration (i.e., progress has stopped), then break
        if len(h_index) == prev_remaining:
            break
        #print('next pair selected')

    # generate summary df
    summary_df = create_summary_df(allele_additions_per_iter, haplotype_additions_per_iter, sample_additions_per_iter)

    # save the summary file
    summary_df.to_csv(os.path.join(output_dir, 'summary_files', 'cross_strat_gRNAs', 'excision_guides','results',gene+'_excision_gRNAs.csv'))
    # note if it finished running
    if summary_df is not None:
        # note if the job completed or failed
        success_fp = os.path.join(
            output_dir,
            "summary_files",
            "cross_strat_gRNAs",
            'excision_guides',
            "logs",
            f"{gene}.DONE_excision"
        )
        with open(success_fp, "w") as f:
            f.write("OK\n")
    else:
        fail_fp = os.path.join(
            output_dir,
            "summary_files",
            "cross_strat_gRNAs",
            'excision_guides',
            "logs",
            f"{gene}.FAIL_excision"
        )
        with open(fail_fp, "w") as f:
            f.write("Failed because summary file is None.\n")


if __name__=='__main__':
    main()
