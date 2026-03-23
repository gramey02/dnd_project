#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
genes_w_guides="$3"
filtered_vcf_dir="$4"
exon_file="$5"
# Local replacement for SGE array index: callers now pass the row number
# explicitly, but we still honor SGE_TASK_ID if this wrapper is reused there.
task_id="${6:-${SGE_TASK_ID:-}}"

if [[ -z "$task_id" ]]; then
    echo "Error: provide a gene row index as argument 6." >&2
    exit 1
fi

# Select the gene assigned to this one local loop iteration.
gene=$(awk -v row="$task_id" 'NR == row {print $1}' "$genes_w_guides")
script="$script_dir/het_combos.py"

python3 "$script" --output_dir "$output_dir" \
    --filtered_vcf_dir "$filtered_vcf_dir" \
    --exon_file "$exon_file" \
    --gene_info "$genes_w_guides" \
    --gene "$gene" \
    --num_samples "$NUM_SAMPLES" \
    --excise_entire_gene "$EXCISE_ENTIRE_GENE"
