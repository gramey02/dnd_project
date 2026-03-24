#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"
execution_utils="$project_root/scripts/utils/execution_mode.sh"

output_dir="$1"
param_file="$2"
source "$param_file"
source "$execution_utils"

# make directory to hold info
mkdir_if_missing() {
  local d="$1"
  if [ -d "$d" ]; then
    echo "Directory '$d' already exists."
  else
    mkdir -p "$d"       # -p creates parents as needed; safe to call repeatedly
    echo "Directory '$d' created."
  fi
}
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/metadata"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/checkpoints"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/results"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/logs"

# merge the information on genes into one file
BASE="$output_dir"
merged_fp="$output_dir/summary_files/cross_strat_gRNAs/metadata/merged_genes_w_valid_guides.txt"

> "$merged_fp"  # truncate/create file

for strat in excision; do
    awk -v s="$strat" '{print $0 "\t" s}' \
        "$BASE/$strat/excavate/het_individuals/metadata/genes_w_valid_guides.txt" \
        >> "$merged_fp"
done

# get the unique genes from the file
unique_genes_file="$output_dir/summary_files/cross_strat_gRNAs/metadata/unique_genes_with_valid_guides_non_excision.txt"
cut -f1 "$merged_fp" | sort -u > "$unique_genes_file"
num_unique_genes=$(wc -l < "$unique_genes_file")
# valid excision pairs
valid_pairs_fp="$output_dir/excision/CommonVars/valid_snp_pairs"

shell_script="$script_dir/excision_guides_array_setup.sh"
run_indexed_jobs "$num_unique_genes" "$shell_script" "$output_dir" "$unique_genes_file" "$param_file" "$valid_pairs_fp"
