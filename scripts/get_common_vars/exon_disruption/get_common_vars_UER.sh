#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../../.." && pwd)"

param_file="$2"
source "$param_file"
exon_file="$4"
output_dir="$1"
total_num_chroms="$3"
# Local replacement for SGE array index: callers now pass the row number
# explicitly, but we still honor SGE_TASK_ID if this wrapper is reused there.
task_id="${5:-${SGE_TASK_ID:-}}"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

# col=$(head -1 $exon_file | tr ',' '\n' | grep -nx "chromosome_name" | cut -d: -f1)
# cut -d',' -f"$col" $exon_file | tail -n +2 | sort -u > $output_dir"/chromosomes/chrom_set.txt"
chrom_set="$resolved_output_base$RUN_NAME/chromosomes/chrom_set.txt"
#chrom_set=$CHROM_SET
if [[ -z "$task_id" ]]; then
  echo "Error: provide a chromosome row index as argument 5." >&2
  exit 1
fi
# Pick the chromosome assigned to this one local loop iteration.
chrom=$(awk -v row="$task_id" 'NR == row {print $1}' "$chrom_set")
echo "$chrom"
# get the file that contains the biallelic snps on the current chromosome
af_file="$AF_FILE_DIR/TGP_chr${chrom}_afs.txt"

# run the indel script
script="$script_dir/get_common_vars.py"
echo "Started identifying common coding vars..."
python3 "$script" --af_limit "$AF_LIMIT" --chrom "$chrom" --af_file "$af_file" --exon_file "$exon_file" --output_dir "$output_dir" --total_num_chroms "$total_num_chroms"
echo "Finished running chromosme ${chrom}."
