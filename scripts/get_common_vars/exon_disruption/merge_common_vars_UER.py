import argparse
import os
import pickle
import pandas as pd

# -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def atomic_pickle_dump(obj, path):
    # write-then-rename so a crash mid-write never leaves a truncated file under the final name
    tmp_path = path + '.tmp'
    with open(tmp_path, 'wb') as f:
        pickle.dump(obj, f)
    os.rename(tmp_path, path)

def atomic_to_csv(df, path, **kwargs):
    tmp_path = path + '.tmp'
    df.to_csv(tmp_path, **kwargs)
    os.rename(tmp_path, path)

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
    parser.add_argument('--output_dir', type=str, required=True, help='Output directory for the pipeline run.')
    parser.add_argument('--total_num_chroms', type=int, required=True, help='Total number of unique chromosomes for gene set of interest.')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    common_vars_savepath = args.output_dir + '/ubiq_region_CommonVars/'
    total_num_chroms = args.total_num_chroms

    summary_files = _collect(common_vars_savepath, 'CommonVars_chr', '_summary.txt', total_num_chroms)
    dict_files = _collect(common_vars_savepath, 'CommonVars_chr', '_dict.pkl', total_num_chroms)

    summary_df = None
    for f in summary_files:
        cur_chrom = f.split("CommonVars_chr")[1].split("_summary.txt")[0]
        cur_df = pd.read_table(common_vars_savepath + f, index_col=0)
        cur_df['chrom'] = cur_chrom
        summary_df = pd.concat([summary_df, cur_df])
    atomic_to_csv(summary_df, common_vars_savepath + 'CommonVars_ALL_summary.txt', sep='\t')
    atomic_to_csv(summary_df, common_vars_savepath + 'CommonVars_ALL_summary_noIDX.txt', sep='\t', header=False, index=False)

    full_dict = {}
    for f in dict_files:
        with open(common_vars_savepath + f, 'rb') as fp:
            full_dict.update(pickle.load(fp))
    atomic_pickle_dump(full_dict, common_vars_savepath + 'CommonVars_ALL_dict.pkl')

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    main()
