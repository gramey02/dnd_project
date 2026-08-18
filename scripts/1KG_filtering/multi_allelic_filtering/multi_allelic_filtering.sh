#!/bin/env bash
#$ -S /bin/bash
#$ -N filter_multiallelic
#$ -cwd
#$ -j y
#$ -o logs/filter_multiallelic_$TASK_ID.log
#$ -l mem_free=4G
#$ -l scratch=20G
#$ -l h_rt=02:00:00
#$ -t 1-23
#$ -tc 23

set -euo pipefail
module load CBI bcftools htslib
IN_DIR="/wynton/home/capra/gramey02/dnd_project/data/biallelic_vcfs_high_coverage"
CSV_DIR="/wynton/home/capra/gramey02/dnd_project/data/multi_allelic_filtering"
OUT_DIR="/wynton/home/capra/gramey02/dnd_project/data/multi_allelic_filtered_vcfs_high_coverage"

mkdir -p "$OUT_DIR" logs

# map SGE_TASK_ID 1-22 -> chrom 1-22, 23 -> X
if [[ "${SGE_TASK_ID}" -eq 23 ]]; then
    chrom="X"
else
    chrom="${SGE_TASK_ID}"
fi

## Set up TMPDIR -> local /scratch, falling back to /tmp if needed
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
id_annotated_vcf="${WORK_DIR}/TGP_chr${chrom}_with_ids.vcf.gz"
out_vcf_local="${WORK_DIR}/TGP_chr${chrom}_biallelicSNPs.vcf.gz"
out_vcf_final="${OUT_DIR}/TGP_chr${chrom}_biallelicSNPs.vcf.gz"

echo "[chr${chrom}] input: $in_vcf"
echo "[chr${chrom}] csv:   $csv"

[[ -f "$in_vcf" ]] || { echo "[chr${chrom}] ERROR: missing $in_vcf" >&2; exit 1; }
[[ -f "$csv" ]]    || { echo "[chr${chrom}] ERROR: missing $csv" >&2; exit 1; }

awk -F',' 'NR>1 && $NF != "" {print $NF}' "$csv" > "$id_list"
n_ids=$(wc -l < "$id_list")
echo "[chr${chrom}] pulled ${n_ids} ids from csv"
[[ "$n_ids" -gt 0 ]] || { echo "[chr${chrom}] WARNING: zero ids parsed" >&2; exit 1; }

# regenerate ID field as CHROM:POS:REF:ALT so it matches the CSV's vcf_id format
# (fixes chrX, which ships with ID="."; no-op on autosomes, which already use this format)
bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' "$in_vcf" -Oz -o "$id_annotated_vcf"
bcftools index -t "$id_annotated_vcf"

bcftools view -i "ID=@${id_list}" "$id_annotated_vcf" -Oz -o "$out_vcf_local"
bcftools index -t "$out_vcf_local"

n_out=$(bcftools index -n "$out_vcf_local")
echo "[chr${chrom}] filtered: ${n_out} variants written locally"

cp "$out_vcf_local" "$out_vcf_final"
cp "${out_vcf_local}.tbi" "${out_vcf_final}.tbi"
echo "[chr${chrom}] done: copied to $out_vcf_final"

[[ -n "${JOB_ID:-}" ]] && qstat -j "$JOB_ID"
