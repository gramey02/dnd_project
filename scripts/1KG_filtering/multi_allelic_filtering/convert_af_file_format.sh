#!/bin/env bash
#$ -S /bin/bash
#$ -N convert_afs
#$ -cwd
#$ -j y
#$ -o logs/convert_afs.log
#$ -l mem_free=2G
#$ -l h_rt=00:30:00

set -euo pipefail

IN_DIR="/wynton/home/capra/gramey02/dnd_project/data/multi_allelic_filtering"
OUT_DIR="/wynton/home/capra/gramey02/dnd_project/data/multiallelic_afs_high_coverage"

mkdir -p "$OUT_DIR" logs

CHROMS=($(seq 1 22) X)

for chrom in "${CHROMS[@]}"; do
    in_csv="${IN_DIR}/chr${chrom}_biallelic_sites_with_rare_multi_0.05.txt"
    out_txt="${OUT_DIR}/TGP_chr${chrom}_afs.txt"

    if [[ ! -f "$in_csv" ]]; then
        echo "[chr${chrom}] WARNING: missing $in_csv, skipping" >&2
        continue
    fi

    # drop header (NR>1), drop index col ($1) and vcf_id col ($NF),
    # keep chr,pos,ref,alt,ac,an,af,afr_af,amr_af,eas_af,eur_af,sas_af
    # output space-delimited, no trailing vcf_id
    awk -F',' 'NR>1 {
        out=$2
        for (i=3; i<NF; i++) out = out" "$i
        print out
    }' "$in_csv" > "$out_txt"

    n_lines=$(wc -l < "$out_txt")
    echo "[chr${chrom}] wrote ${n_lines} lines to $out_txt"
done

echo "Done."
