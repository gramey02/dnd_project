#!/bin/bash
#$ -N find_excision_commonVars
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

output_dir="$1"
param_file="$2"
source "$param_file"
gene_info="$3"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

script="$script_dir/get_common_vars/excision/find_excision_commonVars.py"
exon_file="$4"

python3 "$script" --output_dir "$output_dir" --gene_info "$gene_info" --af_limit "$AF_LIMIT" --af_file_dir "$project_root/$AF_FILE_DIR" --exon_file "$project_root/$exon_file" --nearby_gene_filter "$FILTER_OUT_NEARBY_GENES"
