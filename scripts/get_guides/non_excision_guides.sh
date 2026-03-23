#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

param_file="$project_root/data/params/params.txt"
source "$param_file"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

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
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs"
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/metadata"
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/checkpoints"
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/results"
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/logs"

# merge the information on genes into one file
BASE="$resolved_output_base$RUN_NAME"
merged_fp="$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/metadata/merged_genes_w_valid_guides.txt"

> "$merged_fp"  # truncate/create file

for strat in indels CRISPRoff acceptor_base_edits donor_base_edits; do
    awk -v s="$strat" '{print $0 "\t" s}' \
        "$BASE/$strat/excavate/het_individuals/metadata/genes_w_valid_guides.txt" \
        >> "$merged_fp"
done

# get the unique genes from the file
unique_genes_file="$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/metadata/unique_genes_with_valid_guides_non_excision.txt"
cut -f1 "$merged_fp" | sort -u > "$unique_genes_file"
num_unique_genes=$(wc -l < "$unique_genes_file")

# set up variable to run all strategies together or separately
all_strats_together="False"

# run locally for each gene
shell_script="$script_dir/non_excision_guides_array_setup.sh"
output_dir="$resolved_output_base$RUN_NAME"
for ((task_id=1; task_id<=num_unique_genes; task_id++)); do
    bash "$shell_script" "$output_dir" "$unique_genes_file" "$param_file" "$all_strats_together" "$task_id"
done
