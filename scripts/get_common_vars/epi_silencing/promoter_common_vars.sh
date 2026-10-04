#!/bin/bash
#$ -N promoter_common_vars
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input params
output_dir="$1"
param_file="$2"
source "$param_file"
exon_file="$3"
total_num_chroms="$4"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

chrom_set=$OUTPUT_DIR"/"$RUN_NAME"/chromosomes/chrom_set.txt"
chrom=$(awk -v row=$SGE_TASK_ID 'NR == row {gsub(/"/,"",$1); print $1}' "$chrom_set")

script="$script_dir/get_common_vars/epi_silencing/promoter_common_vars.py"
python3 "$script" --cpg_file "$project_root/$CPG_FILE" \
  --gt_file "$project_root/$exon_file" \
  --intersection "$INTERSECTION" \
  --promoter_ud "$PROMOTER_UD" \
  --promoter_dd "$PROMOTER_DD" \
  --output_dir "$output_dir" \
  --af_limit "$AF_LIMIT" \
  --af_file_dir "$project_root/$AF_FILE_DIR" \
  --gc_threshold "$GC_THRESH" \
  --use_islands "$USE_ISLANDS" \
  --ref_genome_fasta "$project_root/$REF_GENOME_FASTA" \
  --chrom "$chrom" \
  --total_num_chroms "$total_num_chroms"
