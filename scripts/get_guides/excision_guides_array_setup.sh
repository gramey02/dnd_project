#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

output_dir="$1"
unique_genes_file="$2"
param_file="$3"
valid_pairs_fp="$4"
# Local replacement for SGE array index: callers now pass the row number
# explicitly, but we still honor SGE_TASK_ID if this wrapper is reused there.
task_id="${5:-${SGE_TASK_ID:-}}"
source "$param_file"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

if [[ -z "$task_id" ]]; then
    echo "Error: provide a gene row index as argument 5." >&2
    exit 1
fi

# Resolve the one gene assigned to this local loop iteration.
gene=$(awk -v row="$task_id" 'NR == row {print $1}' "$unique_genes_file")
vcf_dir="$resolved_output_base$RUN_NAME/excision/excavate/Guide_filtered_vcfs" # "$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/downsampled_vcfs"

OUTPUT_FILE="$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/results/${gene}_excision_gRNAs.csv"

if [ ! -f "$OUTPUT_FILE" ]; then
    # downsample to the top 258 snp locations with the higest heterozygote frequency
    gRNA_script="$script_dir/excision_guides.py"
    python3 "$gRNA_script" --output_dir "$output_dir" --gene "$gene" --num_samples "$NUM_SAMPLES" --valid_pairs_fp "$valid_pairs_fp" --vcf_dir "$vcf_dir" --max_iter 50
fi
