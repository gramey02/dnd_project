import argparse
import os
import pickle
import pandas as pd
from get_donor_splice_site_vars import atomic_pickle_dump, atomic_to_csv

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
    parser.add_argument('--output_dir', type=str, required=True, help='Location for output files.')
    parser.add_argument('--total_num_chroms', type=int, required=True, help='Total number of unique chromosomes for gene set of interest.')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    output_dir = args.output_dir
    total_num_chroms = args.total_num_chroms

    ubiq_dir = output_dir + '/ubiq_regions/'
    common_vars_savepath = output_dir + '/ubiq_region_CommonVars/'

    # merge universal donor snp regions
    donor_region_files = _collect(ubiq_dir, 'ubiq_donorRegions_chr', '.pkl', total_num_chroms)
    full_donor_regions = {}
    for f in donor_region_files:
        with open(ubiq_dir + f, 'rb') as fp:
            full_donor_regions.update(pickle.load(fp))
    atomic_pickle_dump(full_donor_regions, ubiq_dir + 'ubiq_donorRegions_ALL_chroms.pkl')

    # merge common-var summaries & dicts
    summary_files = _collect(common_vars_savepath, 'CommonVars_chr', '_summary.txt', total_num_chroms)
    dict_files = _collect(common_vars_savepath, 'CommonVars_chr', '_dict.pkl', total_num_chroms)
    summary_df = None
    for f in summary_files:
        cur_df = pd.read_table(common_vars_savepath + f, index_col=0)
        summary_df = pd.concat([summary_df, cur_df])
    atomic_to_csv(summary_df, common_vars_savepath + 'CommonVars_ALL_summary.txt', sep='\t')
    atomic_to_csv(summary_df, common_vars_savepath + 'CommonVars_ALL_summary_noIDX.txt', sep='\t', index=False, header=False)
    full_dict = {}
    for f in dict_files:
        with open(common_vars_savepath + f, 'rb') as fp:
            full_dict.update(pickle.load(fp))
    atomic_pickle_dump(full_dict, common_vars_savepath + 'CommonVars_ALL_dict.pkl')

    # merge base editor summaries
    be_summary_files = _collect(common_vars_savepath, 'base_editor_chr', '_summary.txt', total_num_chroms)
    be_all_df = None
    for f in be_summary_files:
        cur_df = pd.read_table(common_vars_savepath + f)
        be_all_df = pd.concat([be_all_df, cur_df])
    atomic_to_csv(be_all_df, common_vars_savepath + 'base_editor_summary.txt', sep='\t', index=False)

# --------------------------------
if __name__ == '__main__':
    main()
