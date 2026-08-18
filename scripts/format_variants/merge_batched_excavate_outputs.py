import argparse
import os
import re
import shutil
import subprocess
import pandas as pd

METADATA_COLS = [
    'hgnc_symbol', 'chrom', 'lower_coord', 'higher_coord',
    'coords', 'filtering_txt_filepath', 'excavate_input_vcf_filepath'
]

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output_dir', type=str, required=True, help="output directory for the pipeline run.")
    args = parser.parse_args()
    return args

def main():
    args = parse_args()
    output_dir = args.output_dir

    metadata_fp = os.path.join(output_dir, 'excavate', 'input_metadata', 'excavate_batch_metadata.txt')
    if not os.path.exists(metadata_fp) or os.path.getsize(metadata_fp) == 0:
        return

    batch_metadata = pd.read_csv(metadata_fp, sep='\t', header=None, names=METADATA_COLS)
    excavate_outputs_dir = os.path.join(output_dir, 'excavate', 'excavate_outputs')

    # group batch pseudo-genes (e.g. DCC_batch1, DCC_batch2) back under their parent gene (DCC)
    batch_metadata['parent_gene'] = batch_metadata['hgnc_symbol'].apply(lambda g: re.sub(r'_batch\d+$', '', g))

    input_vcfs_dir = os.path.join(output_dir, 'excavate', 'input_vcfs')

    for parent_gene, group in batch_metadata.groupby('parent_gene'):
        batch_output_dirs = [
            os.path.join(excavate_outputs_dir, f'{batch_gene}_output')
            for batch_gene in group['hgnc_symbol']
        ]

        guides_dfs = [
            pd.read_csv(os.path.join(d, 'all_guides.csv'))
            for d in batch_output_dirs if os.path.exists(os.path.join(d, 'all_guides.csv'))
        ]
        summary_dfs = [
            pd.read_csv(os.path.join(d, 'all_guides_summary.csv'))
            for d in batch_output_dirs if os.path.exists(os.path.join(d, 'all_guides_summary.csv'))
        ]

        merged_output_dir = os.path.join(excavate_outputs_dir, f'{parent_gene}_output')
        os.makedirs(merged_output_dir, exist_ok=True)

        if guides_dfs:
            pd.concat(guides_dfs, ignore_index=True).to_csv(
                os.path.join(merged_output_dir, 'all_guides.csv'), index=False
            )
        if summary_dfs:
            pd.concat(summary_dfs, ignore_index=True).to_csv(
                os.path.join(merged_output_dir, 'all_guides_summary.csv'), index=False
            )

        # remove the now-redundant per-batch output dirs so downstream steps
        # that list excavate_outputs/ don't double-count this gene
        for d in batch_output_dirs:
            if os.path.exists(d):
                shutil.rmtree(d)

        # consolidate the per-batch VCFs into one VCF for the parent gene, sorted by position
        batch_vcfs = [
            os.path.join(input_vcfs_dir, f'{batch_gene}_CommonVar_filtered.vcf.gz')
            for batch_gene in group['hgnc_symbol']
        ]
        batch_vcfs = [f for f in batch_vcfs if os.path.exists(f)]
        if batch_vcfs:
            merged_vcf = os.path.join(input_vcfs_dir, f'{parent_gene}_CommonVar_filtered.vcf.gz')
            concat = subprocess.run(
                ['bcftools', 'concat', '-a'] + batch_vcfs,
                stdout=subprocess.PIPE, check=True
            )
            subprocess.run(
                ['bcftools', 'sort', '-Oz', '-o', merged_vcf],
                input=concat.stdout, check=True
            )
            subprocess.run(['bcftools', 'index', '-t', merged_vcf], check=True)

            for f in batch_vcfs:
                os.remove(f)
                tbi = f + '.tbi'
                if os.path.exists(tbi):
                    os.remove(tbi)

if __name__ == '__main__':
    main()
