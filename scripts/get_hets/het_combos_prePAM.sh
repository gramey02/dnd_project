#!/bin/bash
#$ -N het_combos
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
gene_info="$3"
filtered_vcf_dir="$4"
valid_pairs_fp="$5"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# set up directories and variables
mkdir -p $output_dir # make output directory if it doesn't already exist
gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $gene_info)

# set script and run
script="$script_dir/get_hets/het_combos_prePAM.py"
python3 "$script" \
    --output_dir "$output_dir" \
    --filtered_vcf_dir "$filtered_vcf_dir" \
    --num_samples "$NUM_SAMPLES" \
    --valid_pairs_fp "$valid_pairs_fp" \
    --gene $gene

