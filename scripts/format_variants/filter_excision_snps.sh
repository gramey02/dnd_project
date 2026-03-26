#!/bin/bash
#$ -N filter_excision_snps
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
cv_dict_filepath="$3"
exon_file="$4"
common_var_genes="$5"

script="$script_dir/filter_excision_snps.py"

gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $common_var_genes)

python3 $script --output_dir $output_dir \
    --gene $gene \
    --cv_dict_filepath $cv_dict_filepath \
    --exon_file $exon_file
