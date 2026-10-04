#!/bin/bash
#$ -N merge_common_vars_UER
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input params
output_dir="$1"
param_file="$2"
source "$param_file"
total_num_chroms="$3"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

script="$script_dir/get_common_vars/exon_disruption/merge_common_vars_UER.py"
python3 "$script" --output_dir "$output_dir" --total_num_chroms "$total_num_chroms"
