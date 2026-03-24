#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

# excavate runs very quickly, so can loop through genes for now, but in the future could also be set up as an array job

# Activate toolchains only when the current environment supports them.
if command -v module >/dev/null 2>&1; then
    module load CBI miniforge3 bcftools
fi
if command -v conda >/dev/null 2>&1; then
    eval "$(conda shell.bash hook)"
    conda activate excavate
fi

# load input vars
output_dir=$1
param_file=$2
source $param_file
gene_info=$3
# Local replacement for SGE array index: callers now pass the row number
# explicitly, but we still honor SGE_TASK_ID if this wrapper is reused there.
task_id="${4:-${SGE_TASK_ID:-}}"

# set variables for inputs and outputs
if [[ "$REF_GENOME_FASTA" = /* ]]; then
    ref_genome_fasta="$REF_GENOME_FASTA"
else
    ref_genome_fasta="$project_root/$REF_GENOME_FASTA"
fi
if [[ -z "$task_id" ]]; then
    echo "Error: provide a gene row index as argument 4." >&2
    exit 1
fi

# Resolve the one gene/chromosome record assigned to this local iteration.
gene=$(awk -v row="$task_id" 'NR == row {print $1}' "$gene_info")
chrom=$(awk -v row="$task_id" 'NR == row {print $2}' "$gene_info")
gene_coords=$(awk -v row="$task_id" 'NR == row {print $5}' "$gene_info")
common_var_vcf=$output_dir"/excavate/input_vcfs/${gene}_CommonVar_filtered.vcf.gz"
if [[ "$CHROM_FASTA_FILEPATH" = /* ]]; then
    chrom_fasta="${CHROM_FASTA_FILEPATH}hg38_${chrom}.fa"
else
    chrom_fasta="$project_root/${CHROM_FASTA_FILEPATH}hg38_${chrom}.fa"
fi

# make output file directory below
output_files=$output_dir"/excavate/excavate_outputs/${gene}_output"
mkdir -p "$output_files"


# run excavate for each gene below
if [[ "$EXCAVATE_SCRIPT" = /* ]]; then
    excavate_script="$EXCAVATE_SCRIPT"
else
    excavate_script="$project_root/$EXCAVATE_SCRIPT"
fi
python3 "$excavate_script" generate "$common_var_vcf" population "$chrom_fasta" "$ref_genome_fasta" "$gene_coords" --cas SpCas9 --summary -o "$output_files"
