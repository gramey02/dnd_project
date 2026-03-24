#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# input vars
output_dir="$1"
unique_genes_file="$2"
param_file="$3"
source "$param_file"
all_strats_together="$4"
# Local replacement for SGE array index: callers now pass the row number
# explicitly, but we still honor SGE_TASK_ID if this wrapper is reused there.
task_id="${5:-${SGE_TASK_ID:-}}"

# Resolve the one gene assigned to this local loop iteration.
if [[ -z "$task_id" ]]; then
    echo "Error: provide a gene row index as argument 5." >&2
    exit 1
fi

gene=$(awk -v row="$task_id" 'NR == row {print $1}' "$unique_genes_file")

gRNA_script="$script_dir/non_excision_guides.py"
python3 "$gRNA_script" --output_dir "$output_dir" --gene "$gene" --num_samples "$NUM_SAMPLES" --all_strats_together "$all_strats_together"
