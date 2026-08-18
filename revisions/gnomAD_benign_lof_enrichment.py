import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import pickle
import random
import subprocess, os
import hail as hl
import argparse

#parse the arguments passed into the script
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--dnd_file', type = str, required = True, help = 'Filtered transcript file.')
    parser.add_argument('--gene_index', type=int, required=True, help='Row index into gene_spans (0-based).')
    parser.add_argument('--ht_path', type=str, required=True, help='Path to pre-pruned local Hail table.')
    parser.add_argument('--output_dir', type=str, required=True, help='Output file path.')
    args = parser.parse_args()
    return args

def is_target_tc(tc, gene):
    return (tc.is_canonical & (tc.gene_symbol == gene) & (tc.lof == 'HC'))

def clinsig_terms(cv_array):
    return hl.or_else(
        hl.flatmap(lambda cv: cv.clin_sig, cv_array),
        hl.empty_array(hl.tstr)
    ).map(lambda s: s.lower())

CATEGORIES = ['pathogenic', 'likely_pathogenic', 'uncertain_significance',
              'likely_benign', 'benign', 'association', 'not_provided']
ALL_CATEGORIES = CATEGORIES + ['other_or_unclassified']


def main():
    args = parse_args()

    dnd = pd.read_csv(args.dnd_file, dtype={'chromosome_name': 'str'})
    chrX_genes = set(dnd[dnd['chromosome_name'] == 'X']['hgnc_symbol'].unique())

    gene_spans = dnd.groupby('hgnc_symbol').agg(
        chrom=('chromosome_name', 'first'),
        start=('start_position', 'min'),
        end=('end_position', 'max')
    ).reset_index()

    # pick this task's single gene by index
    r = gene_spans.iloc[args.gene_index]
    gene = r.hgnc_symbol

    # skip re-running if output already exists (useful for resubmitted/retried tasks)
    out_path = os.path.join(args.output_dir, f"{gene}.csv")
    if os.path.exists(out_path):
        print(f"Output already exists for {gene}, skipping.")
        return

    hl.init(log=f"/tmp/hail_{gene}.log")
    hl.default_reference("GRCh38")

    # read the pre-pruned local table (no S3/internet needed on compute nodes)
    ht = hl.read_table(args.ht_path)

    chrom_str = str(r.chrom)
    if not chrom_str.startswith('chr'):
        chrom_str = 'chr' + chrom_str
    interval = hl.parse_locus_interval(f"{chrom_str}:{r.start}-{r.end}", reference_genome='GRCh38')
    ht_region = hl.filter_intervals(ht, [interval])

    ht_lof = ht_region.filter(
        hl.any(lambda tc: is_target_tc(tc, gene), ht_region.transcript_consequences)
    )

    ht_lof = ht_lof.annotate(
        ac_total = hl.or_else(ht_lof.exome.freq.all.ac, 0) + hl.or_else(ht_lof.genome.freq.all.ac, 0),
        hom_total = hl.or_else(ht_lof.exome.freq.all.homozygote_count, hl.int64(0)) +
                    hl.or_else(ht_lof.genome.freq.all.homozygote_count, hl.int64(0)),
        hemi_total = hl.or_else(ht_lof.exome.freq.all.hemizygote_count, 0) +
                     hl.or_else(ht_lof.genome.freq.all.hemizygote_count, 0),
    )
    ht_lof = ht_lof.annotate(het_total = ht_lof.ac_total - (2 * ht_lof.hom_total) - ht_lof.hemi_total)

    is_chrX_gene = gene in chrX_genes
    ht_het = ht_lof.filter(ht_lof.hemi_total > 0) if is_chrX_gene else ht_lof.filter(ht_lof.het_total > 0)

    ht_het = ht_het.annotate(
        clinsig_terms = clinsig_terms(ht_het.exome.vep115.colocated_variants).extend(
                        clinsig_terms(ht_het.genome.vep115.colocated_variants))
    )
    ht_het = ht_het.annotate(
        clinsig_category = (
            hl.case()
            .when(ht_het.clinsig_terms.any(lambda s: s.contains('pathogenic') & ~s.contains('likely')), 'pathogenic')
            .when(ht_het.clinsig_terms.any(lambda s: s.contains('likely_pathogenic')), 'likely_pathogenic')
            .when(ht_het.clinsig_terms.any(lambda s: s.contains('uncertain')), 'uncertain_significance')
            .when(ht_het.clinsig_terms.any(lambda s: s.contains('likely_benign')), 'likely_benign')
            .when(ht_het.clinsig_terms.any(lambda s: s.contains('benign') & ~s.contains('likely')), 'benign')
            .when(ht_het.clinsig_terms.any(lambda s: s.contains('association')), 'association')
            .when(ht_het.clinsig_terms.any(lambda s: s.contains('not_provided')), 'not_provided')
            .default('other_or_unclassified')
        )
    )

    gene_counts_df = ht_het.group_by(ht_het.clinsig_category).aggregate(n=hl.agg.count()).to_pandas()

    gene_counts_dict = {cat: 0 for cat in ALL_CATEGORIES}
    gene_counts_dict.update(dict(zip(gene_counts_df['clinsig_category'], gene_counts_df['n'])))

    out_df = pd.DataFrame([gene_counts_dict], index=[gene])
    out_df.index.name = 'gene_symbol'
    os.makedirs(args.output_dir, exist_ok=True)
    out_df.to_csv(out_path)

    hl.stop()

if __name__ == '__main__':
    main()