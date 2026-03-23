#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

# Load cluster modules only if the current environment supports `module`.
if command -v module >/dev/null 2>&1; then
    module load CBI bcftools
fi

# load input arguments
output_dir=$1
param_file=$2
source $param_file
gene_info=$3
# Local replacement for SGE array index: callers now pass the row number
# explicitly, but we still honor SGE_TASK_ID if this wrapper is reused there.
task_id="${4:-${SGE_TASK_ID:-}}"


# separate out inputs and outputs
if [[ -z "$task_id" ]]; then
    echo "Error: provide a gene row index as argument 4." >&2
    exit 1
fi

# Select the gene/chromosome pair assigned to this one local loop iteration.
chrom=$(awk -v row="$task_id" 'NR == row {print $2}' "$gene_info")
gene=$(awk -v row="$task_id" 'NR == row {print $1}' "$gene_info")
positions=$output_dir"/excavate/Guide_locs/${gene}_Guide_locs.txt"
input_vcf="$PHASED_1000G_VCF_DIR/ALL.chr${chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
filtered_vcf=$output_dir"/excavate/Guide_filtered_vcfs/${gene}_guide_filtered.vcf.gz"

# filter with bcftools
bcftools view -T "$positions" "$input_vcf" -Oz -o "$filtered_vcf"
# unzip so you can run analysis on the file later
gzip -dk "$filtered_vcf"









# non-parallelized
# # test:
# gene_info="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/indel_PAM_filtered_vars/dhs_small.txt"
# # full:
# #gene_info="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/indel_PAM_filtered_vars/dhs_indel_PAM_2025_04_22_noIDX.txt"

# edit_strategy="indel"

# while read gene chrom start_pos end_pos gene_coords; do
#     echo "Started ${gene} filtering..."
#     # assign variables for the current gene
#     positions="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/${edit_strategy}_PAM_filtered_vars/${gene}_${edit_strategy}_PAM_positions.txt"
#     input_vcf="/wynton/group/capra/data/wynton_databases/1000_genomes/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
#     filtered_vcf="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/${edit_strategy}_PAM_filtered_vars/${gene}.vcf.gz"
#     bcftools view -T $positions $input_vcf -Oz -o $filtered_vcf
#     gzip -dk $filtered_vcf
#     echo "Finished ${gene} filtering."
# done < $gene_info
