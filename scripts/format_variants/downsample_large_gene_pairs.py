import argparse
import os
import pickle
import pandas as pd

VCF_FIXED_COLS = ['chr', 'pos', 'rsid', 'ref', 'alt', 'qual', 'filter', 'info', 'format']

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genes_w_guides_large', type=str, required=True, help='File listing genes routed to the large het_combos job.')
    parser.add_argument('--valid_pairs_fp', type=str, required=True, help='Directory of original <gene>_valid_snp_pairs.pkl files.')
    parser.add_argument('--filtered_vcf_dir', type=str, required=True, help='Directory of guide-filtered vcfs.')
    parser.add_argument('--num_samples', type=int, required=True, help='Number of samples in the population.')
    parser.add_argument('--ds_threshold', type=int, required=True, help='Max number of SNP pairs to keep per gene.')
    parser.add_argument('--output_dir', type=str, required=True, help='Directory to write downsampled <gene>_valid_snp_pairs.pkl files to.')
    args = parser.parse_args()
    return args

def get_snp_het_counts(vcf_fp, num_samples):
    """Return a dict mapping each SNP position to its heterozygote count."""
    sample_cols = ["sample" + str(s) for s in range(1, num_samples + 1)]
    vcf = pd.read_table(vcf_fp, comment='#', header=None)
    vcf.columns = VCF_FIXED_COLS + sample_cols

    het_mask = vcf[sample_cols].isin(["0|1", "1|0"])
    het_counts = het_mask.sum(axis=1)
    return dict(zip(vcf['pos'], het_counts))

def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    with open(args.genes_w_guides_large) as f:
        genes = [line.split('\t')[0].strip() for line in f if line.strip()]

    for gene in genes:
        pairs_fp = os.path.join(args.valid_pairs_fp, f'{gene}_valid_snp_pairs.pkl')
        if not os.path.exists(pairs_fp):
            continue

        with open(pairs_fp, 'rb') as fp:
            pairs = list(pickle.load(fp))

        out_fp = os.path.join(args.output_dir, f'{gene}_valid_snp_pairs.pkl')

        if len(pairs) <= args.ds_threshold:
            with open(out_fp, 'wb') as fp:
                pickle.dump(pairs, fp)
            continue

        vcf_fp = os.path.join(args.filtered_vcf_dir, f'{gene}_guide_filtered.vcf')
        het_counts = get_snp_het_counts(vcf_fp, args.num_samples)

        # score each pair by the combined heterozygote count of its two snps, so pairs
        # more likely to matter for the greedy algorithm are preferentially kept
        scored_pairs = sorted(
            pairs,
            key=lambda pair: het_counts.get(pair[0], 0) + het_counts.get(pair[1], 0),
            reverse=True,
        )
        downsampled_pairs = scored_pairs[:args.ds_threshold]

        with open(out_fp, 'wb') as fp:
            pickle.dump(downsampled_pairs, fp)

if __name__ == '__main__':
    main()
