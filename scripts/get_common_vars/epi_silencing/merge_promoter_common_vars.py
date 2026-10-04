import argparse
import os
import pickle
import pandas as pd
from promoter_common_vars import atomic_pickle_dump, atomic_to_csv

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def _collect(directory, prefix, suffix, total_num_chroms):
    files = [f for f in os.listdir(directory) if f.startswith(prefix) and f.endswith(suffix)]
    if len(files) != total_num_chroms:
        raise RuntimeError(
            f"Expected {total_num_chroms} files matching '{prefix}*{suffix}' in {directory}, found {len(files)}: {files}"
        )
    return files

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory to save information to.')
    parser.add_argument('--total_num_chroms', type=int, required=True, help='Total number of unique chromosomes for gene set of interest.')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    output_dir = args.output_dir
    total_num_chroms = args.total_num_chroms

    ubiq_dir = output_dir + '/ubiq_regions/'
    common_vars_dir = output_dir + '/ubiq_region_CommonVars/'

    # merge CpG island overlap dicts
    cpg_files = _collect(ubiq_dir, 'dhs_CpGIsland_Overlap_chr', '_dict.pkl', total_num_chroms)
    full_cpg_overlap_dict = {}
    for f in cpg_files:
        with open(ubiq_dir + f, 'rb') as fp:
            full_cpg_overlap_dict.update(pickle.load(fp))
    atomic_pickle_dump(full_cpg_overlap_dict, ubiq_dir + 'dhs_CpGIsland_Overlap_ALL_dict.pkl')

    # merge shared promoter region dicts
    promoter_files = _collect(ubiq_dir, 'ubiq_promoters_chr', '.pkl', total_num_chroms)
    full_shared_promoter_regions = {}
    for f in promoter_files:
        with open(ubiq_dir + f, 'rb') as fp:
            full_shared_promoter_regions.update(pickle.load(fp))
    atomic_pickle_dump(full_shared_promoter_regions, ubiq_dir + 'ubiq_promoters_ALL_chroms.pkl')

    # merge common-var summaries & dicts
    summary_files = _collect(common_vars_dir, 'CommonVars_chr', '_summary.txt', total_num_chroms)
    dict_files = _collect(common_vars_dir, 'CommonVars_chr', '_dict.pkl', total_num_chroms)
    summary_df = None
    for f in summary_files:
        cur_df = pd.read_table(common_vars_dir + f, index_col=0)
        summary_df = pd.concat([summary_df, cur_df])
    atomic_to_csv(summary_df, common_vars_dir + 'CommonVars_ALL_summary.txt', sep='\t')
    summary_df_small = summary_df[['hgnc_symbol', 'num_common_vars_in_shared_promoters', 'chromosome_name']]
    atomic_to_csv(summary_df_small, common_vars_dir + 'CommonVars_ALL_summary_noIDX.txt', sep='\t', header=None, index=None)
    full_dict = {}
    for f in dict_files:
        with open(common_vars_dir + f, 'rb') as fp:
            full_dict.update(pickle.load(fp))
    atomic_pickle_dump(full_dict, common_vars_dir + 'CommonVars_ALL_dict.pkl')

# -----------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
