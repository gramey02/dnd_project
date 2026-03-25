#!/bin/bash
#$ -N get_common_vars_UER
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ../../../logs/out/get_common_vars_UER.out
#$ -e ../../../logs/err/get_common_vars_UER.err

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../../.." && pwd)"

param_file="$2"
source "$param_file"
exon_file="$4"
output_dir="$1"
total_num_chroms="$3"

# col=$(head -1 $exon_file | tr ',' '\n' | grep -nx "chromosome_name" | cut -d: -f1)
# cut -d',' -f"$col" $exon_file | tail -n +2 | sort -u > $output_dir"/chromosomes/chrom_set.txt"
chrom_set=$OUTPUT_DIR"/"$RUN_NAME"/chromosomes/chrom_set.txt"
#chrom_set=$CHROM_SET
chrom=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $chrom_set)

# run the indel script
script="$script_dir/get_common_vars.py"
python3 "$script" --af_limit "$AF_LIMIT" --chrom "$chrom" --af_file "$af_file" --exon_file "$exon_file" --output_dir "$output_dir" --total_num_chroms "$total_num_chroms"
