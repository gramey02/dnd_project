#!/bin/bash
#$ -N get_targeted_hets
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input arguments
output_dir="$1"
param_file="$2"
source "$param_file"
gene_info="$3"
excavate_output_dir="$4"
filtered_vcf_dir="$5"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# run het individual calculation script
# (script will save both the number of hets and the identifiers of the hets
# to be compared to other editing methods)
script="$script_dir/get_hets/get_targeted_hets.py"
python3 "$script" --output_dir "$output_dir" \
    --gene_info "$gene_info" \
    --excavate_output_dir "$excavate_output_dir" \
    --filtered_vcf_dir "$filtered_vcf_dir"
