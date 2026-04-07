#!/bin/bash
#$ -N get_donor_splice_site_vars
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
exon_file="$3"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

script="$script_dir/get_common_vars/ss_disruption/get_donor_splice_site_vars.py"

# Run the script below
python3 "$script" --exon_file "$project_root/$exon_file" \
  --af_limit "$AF_LIMIT" \
  --af_file_dir "$project_root/$AF_FILE_DIR" \
  --editing_window_size "$EDITING_WINDOW_SIZE" \
  --donor_snp_region "$DONOR_SNP_REGION" \
  --output_dir "$output_dir"
