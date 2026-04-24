#!/bin/bash
#$ -N non_excision_guides_array_setup
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input args
output_dir="$1"
unique_genes_file="$2"
param_file="$3"
source "$param_file"
all_strats_together="$4"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# run greedy algorithm on one gene per array job
gene=$(awk -v row="$SGE_TASK_ID" 'NR == row {print $1}' "$unique_genes_file")
gRNA_script="$script_dir/get_guides/non_excision_guides.py"
python3 "$gRNA_script" --output_dir "$output_dir" --gene "$gene" --num_samples "$NUM_SAMPLES" --all_strats_together "$all_strats_together"
