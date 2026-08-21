#!/bin/env bash
#$ -S /bin/bash
#$ -N validate_multiallelic
#$ -cwd
#$ -j y
#$ -o logs/validate_multiallelic.log
#$ -l mem_free=2G
#$ -l scratch=20G
#$ -l h_rt=00:30:00

set -uo pipefail   # no -e: keep checking all chroms even if one fails
module load CBI bcftools htslib
CSV_DIR="/wynton/home/capra/gramey02/dnd_project/data/multi_allelic_filtering"
OUT_DIR="/wynton/home/capra/gramey02/dnd_project/data/multi_allelic_filtered_vcfs"

mkdir -p logs

## 0. Set up TMPDIR -> local /scratch, falling back to /tmp if needed
if [[ -z "${TMPDIR:-}" ]]; then
  if [[ -d /scratch ]]; then TMPDIR="/scratch/$USER"; else TMPDIR="/tmp/$USER"; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

WORK_DIR="${TMPDIR}/validate_multiallelic"
mkdir -p "$WORK_DIR"

CHROMS=($(seq 1 22) X)

echo "chrom,csv_rows,vcf_rows,match"
echo "-------------------------------"

n_mismatch=0
n_missing=0

for chrom in "${CHROMS[@]}"; do
    echo "[chr${chrom}] checking..." >&2

    csv="${CSV_DIR}/chr${chrom}_biallelic_sites_with_rare_multi_0.05.txt"
    vcf="${OUT_DIR}/TGP_chr${chrom}_biallelicSNPs_filtered.vcf.gz"

    if [[ ! -f "$csv" ]]; then
        echo "chr${chrom},MISSING_CSV,-,FAIL"
        n_missing=$((n_missing+1))
        continue
    fi
    if [[ ! -f "$vcf" ]]; then
        echo "chr${chrom},-,MISSING_VCF,FAIL"
        n_missing=$((n_missing+1))
        continue
    fi

    csv_ids="${WORK_DIR}/chr${chrom}_csv_ids.txt"
    awk -F',' 'NR>1 && $NF != "" {print $NF}' "$csv" | LC_ALL=C sort -T "$WORK_DIR" -u > "$csv_ids"
    csv_rows=$(wc -l < "$csv_ids")

    if [[ -f "${vcf}.tbi" || -f "${vcf}.csi" ]]; then
        vcf_rows=$(bcftools index -n "$vcf")
    else
        vcf_rows=$(zcat "$vcf" | grep -vc '^#')
    fi

    if [[ "$csv_rows" -eq "$vcf_rows" ]]; then
        match="OK"
    else
        match="MISMATCH"
        n_mismatch=$((n_mismatch+1))
    fi

    echo "chr${chrom},${csv_rows},${vcf_rows},${match}"
done

echo "-------------------------------"
if [[ "$n_mismatch" -eq 0 && "$n_missing" -eq 0 ]]; then
    echo "All chromosomes matched exactly."
else
    echo "${n_mismatch} chromosome(s) with row count mismatch, ${n_missing} missing file(s). See above."
fi

[[ -n "${JOB_ID:-}" ]] && qstat -j "$JOB_ID"
