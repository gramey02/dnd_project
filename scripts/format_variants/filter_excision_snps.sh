#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
cv_dict_filepath="$3"
exon_file="$4"
common_var_genes="$5"
task_id="${6:-${SGE_TASK_ID:-}}"

script="$script_dir/filter_excision_snps.py"

if [[ -z "$task_id" ]]; then
    echo "Error: provide a gene row index as argument 6." >&2
    exit 1
fi

gene=$(awk -v row="$task_id" 'NR == row {print $1}' "$common_var_genes")

python3 "$script" --output_dir "$output_dir" \
    --gene "$gene" \
    --cv_dict_filepath "$cv_dict_filepath" \
    --exon_file "$exon_file"
