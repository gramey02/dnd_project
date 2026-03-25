#!/bin/bash
#$ -N generate_filtered_vcfs
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ../../logs/out/generate_filtered_vcfs.out
#$ -e ../../logs/err/generate_filtered_vcfs.err

set -euo pipefail

module load CBI bcftools

# file that contains metadata on genes and coordinates
output_dir="$1"
param_file="$2"
source "$param_file"
gene_info="$3"
task_id="${4:-${SGE_TASK_ID:-}}"

if [[ -z "$task_id" ]]; then
    echo "Error: provide a gene row index as argument 4 or run as an SGE array job." >&2
    exit 1
fi

# set up to run an array job
cur_gene=$(awk -v row="$task_id" 'NR == row {print $1}' "$gene_info")
cur_chrom=$(awk -v row="$task_id" 'NR == row {print $2}' "$gene_info")
# get current gene's common variant position file
common_vars_pos="$output_dir/excavate/CommonVar_locs/${cur_gene}_CommonVar_locs.txt"

# output directory for the current analysis
#OUTPUT=
# filter the vcf on the current chromosome accordingly to these variants
output_vcf="$output_dir/excavate/input_vcfs/${cur_gene}_CommonVar_filtered.vcf.gz"
# assign biallelic snp file based on current chromosome
biallelic_snps="$BIALLELIC_SNPS_DIR/TGP_chr${cur_chrom}_biallelicSNPs.vcf.gz"

# run filtering
bcftools view -T "$common_vars_pos" "$biallelic_snps" -Oz -o "$output_vcf"
# create index file
bcftools index -t "$output_vcf"
