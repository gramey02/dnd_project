import argparse
import os
import pickle

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--genes_w_guides', type=str, required=True)
    parser.add_argument('--valid_pairs_fp', type=str, required=True)
    parser.add_argument('--threshold', type=int, required=True)
    parser.add_argument('--normal_out', type=str, required=True)
    parser.add_argument('--large_out', type=str, required=True)
    return parser.parse_args()

def main():
    args = parse_args()
    normal_rows, large_rows = [], []

    with open(args.genes_w_guides) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            gene = line.split('\t')[0]
            pairs_fp = os.path.join(args.valid_pairs_fp, f'{gene}_valid_snp_pairs.pkl')
            n_pairs = 0
            if os.path.exists(pairs_fp):
                with open(pairs_fp, 'rb') as fp:
                    n_pairs = len(pickle.load(fp))
            (large_rows if n_pairs > args.threshold else normal_rows).append(line)

    for rows, out_fp in [(normal_rows, args.normal_out), (large_rows, args.large_out)]:
        with open(out_fp, 'w') as f:
            f.write('\n'.join(rows) + ('\n' if rows else ''))

if __name__ == '__main__':
    main()
