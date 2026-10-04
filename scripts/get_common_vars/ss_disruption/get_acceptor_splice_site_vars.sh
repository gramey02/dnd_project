#!/bin/bash
#$ -N get_acceptor_splice_site_vars
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
exon_file="$3"
total_num_chroms="$4"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

chrom_set=$OUTPUT_DIR"/"$RUN_NAME"/chromosomes/chrom_set.txt"
chrom=$(awk -v row=$SGE_TASK_ID 'NR == row {gsub(/"/,"",$1); print $1}' "$chrom_set")

script="$script_dir/get_common_vars/ss_disruption/get_acceptor_splice_site_vars.py"

# Run the script below
python3 "$script" --exon_file "$project_root/$exon_file" \
  --af_limit "$AF_LIMIT" \
  --af_file_dir "$project_root/$AF_FILE_DIR" \
  --editing_window_size "$EDITING_WINDOW_SIZE" \
  --acceptor_snp_region "$ACCEPTOR_SNP_REGION" \
  --output_dir "$output_dir" \
  --base_editor "$BASE_EDITOR" \
  --chrom "$chrom" \
  --total_num_chroms "$total_num_chroms"
