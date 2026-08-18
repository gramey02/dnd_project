import argparse
import os
import pandas as pd

METADATA_COLS = [
    'hgnc_symbol', 'chrom', 'lower_coord', 'higher_coord',
    'coords', 'filtering_txt_filepath', 'excavate_input_vcf_filepath'
]

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type=str, required=True, help="output directory for the pipeline run.")
    parser.add_argument('--ds_threshold', type=int, required=True, help='max number of variants to run through EXCAVATE in a single batch')
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    output_dir = args.output_dir
    batch_threshold = args.ds_threshold

    common_var_locs_dir = os.path.join(output_dir, 'excavate', 'CommonVar_locs')
    metadata_dir = os.path.join(output_dir, 'excavate', 'input_metadata')
    metadata_fp = os.path.join(metadata_dir, 'excavate_run_metadata.txt')

    metadata = pd.read_csv(metadata_fp, sep='\t', header=None, names=METADATA_COLS)

    unbatched_rows = []
    batch_rows = []

    for _, gene_row in metadata.iterrows():
        gene = gene_row['hgnc_symbol']
        locs_fp = os.path.join(common_var_locs_dir, f'{gene}_CommonVar_locs.txt')
        locs_df = pd.read_csv(locs_fp, sep='\t', header=None, names=['chrom', 'pos'])

        if len(locs_df) <= batch_threshold:
            unbatched_rows.append(gene_row)
            continue

        # split this gene's variants into sequential batches of at most batch_threshold variants each
        for batch_num, start in enumerate(range(0, len(locs_df), batch_threshold), start=1):
            batch_locs_df = locs_df.iloc[start:start + batch_threshold]
            batch_gene = f'{gene}_batch{batch_num}'

            batch_locs_fp = os.path.join(common_var_locs_dir, f'{batch_gene}_CommonVar_locs.txt')
            batch_locs_df.to_csv(batch_locs_fp, sep='\t', header=False, index=False)

            batch_row = gene_row.copy()
            batch_row['hgnc_symbol'] = batch_gene
            batch_row['filtering_txt_filepath'] = batch_locs_fp
            batch_row['excavate_input_vcf_filepath'] = os.path.join(
                output_dir, 'excavate', 'input_vcfs', f'{batch_gene}_CommonVar_filtered.vcf.gz'
            )
            batch_rows.append(batch_row)

    unbatched_metadata = pd.DataFrame(unbatched_rows, columns=METADATA_COLS)
    unbatched_metadata.to_csv(
        os.path.join(metadata_dir, 'excavate_run_metadata_unbatched.txt'),
        sep='\t', header=False, index=False
    )

    batch_metadata = pd.DataFrame(batch_rows, columns=METADATA_COLS)
    batch_metadata.to_csv(
        os.path.join(metadata_dir, 'excavate_batch_metadata.txt'),
        sep='\t', header=False, index=False
    )

if __name__ == '__main__':
    main()
