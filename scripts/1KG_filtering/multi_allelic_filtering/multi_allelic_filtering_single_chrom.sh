#!/bin/env bash
#$ -S /bin/bash
#$ -N filter_multiallelic2
#$ -cwd
#$ -j y
#$ -o logs/filter_multiallelic2.log
#$ -l mem_free=4G
#$ -l scratch=20G
#$ -l h_rt=02:00:00

set -euo pipefail
module load CBI bcftools htslib
# ---- usage ----
if [[ $# -ne 1 ]]; then
    echo "Usage: qsub filter_multiallelic_single.sge.sh <chrom>" >&2
    echo "  e.g.: qsub filter_multiallelic_single.sge.sh 1" >&2
    echo "        qsub filter_multiallelic_single.sge.sh X" >&2
    exit 1
fi

chrom="$1"

# ---- config ----
IN_DIR="/wynton/home/capra/gramey02/dnd_project/data/biallelic_vcfs_high_coverage"
CSV_DIR="/wynton/home/capra/gramey02/dnd_project/data/multi_allelic_filtering"
OUT_DIR="/wynton/home/capra/gramey02/dnd_project/data/multi_allelic_filtered_vcfs_high_coverage"

mkdir -p "$OUT_DIR" logs

## 0. Set up TMPDIR -> local /scratch, falling back to /tmp if needed
if [[ -z "${TMPDIR:-}" ]]; then
  if [[ -d /scratch ]]; then TMPDIR="/scratch/$USER"; else TMPDIR="/tmp/$USER"; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

WORK_DIR="${TMPDIR}/filter_multiallelic_chr${chrom}"
mkdir -p "$WORK_DIR"

in_vcf="${IN_DIR}/TGP_chr${chrom}_biallelicSNPs.vcf.gz"
csv="${CSV_DIR}/chr${chrom}_biallelic_sites_with_rare_multi_0.05.txt"
id_list="${WORK_DIR}/chr${chrom}_ids.txt"
out_vcf_local="${WORK_DIR}/TGP_chr${chrom}_biallelicSNPs.vcf.gz"
out_vcf_final="${OUT_DIR}/TGP_chr${chrom}_biallelicSNPs.vcf.gz"

echo "[chr${chrom}] TMPDIR:        $TMPDIR"
echo "[chr${chrom}] input:         $in_vcf"
echo "[chr${chrom}] csv:           $csv"
echo "[chr${chrom}] local output:  $out_vcf_local"
echo "[chr${chrom}] final output:  $out_vcf_final"

# sanity checks
[[ -f "$in_vcf" ]] || { echo "[chr${chrom}] ERROR: missing $in_vcf" >&2; exit 1; }
[[ -f "$csv" ]]    || { echo "[chr${chrom}] ERROR: missing $csv" >&2; exit 1; }

# extract vcf_id (last column) from csv, skipping header, dropping any blank lines
awk -F',' 'NR>1 && $NF != "" {print $NF}' "$csv" > "$id_list"

n_ids=$(wc -l < "$id_list")
echo "[chr${chrom}] pulled ${n_ids} ids from csv"

if [[ "$n_ids" -eq 0 ]]; then
    echo "[chr${chrom}] WARNING: zero ids parsed, skipping bcftools call" >&2
    exit 1
fi

# filter vcf to only these ids -- write to local scratch first
bcftools view -i "ID=@${id_list}" "$in_vcf" -Oz -o "$out_vcf_local"
bcftools index -t "$out_vcf_local"

n_out=$(bcftools index -n "$out_vcf_local")
echo "[chr${chrom}] filtered: ${n_out} variants written locally"

# copy final files back to persistent storage
cp "$out_vcf_local" "$out_vcf_final"
cp "${out_vcf_local}.tbi" "${out_vcf_final}.tbi"

echo "[chr${chrom}] done: copied to $out_vcf_final"

# ---- end-of-job summary ----
[[ -n "${JOB_ID:-}" ]] && qstat -j "$JOB_ID"
