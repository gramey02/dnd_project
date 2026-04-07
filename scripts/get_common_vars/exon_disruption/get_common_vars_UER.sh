#!/bin/bash
#$ -N get_common_vars_UER
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

param_file="$2"
source "$param_file"
exon_file="$4"
output_dir="$1"
total_num_chroms="$3"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# col=$(head -1 $exon_file | tr ',' '\n' | grep -nx "chromosome_name" | cut -d: -f1)
# cut -d',' -f"$col" $exon_file | tail -n +2 | sort -u > $output_dir"/chromosomes/chrom_set.txt"
chrom_set=$OUTPUT_DIR"/"$RUN_NAME"/chromosomes/chrom_set.txt"
#chrom_set=$CHROM_SET
chrom=$(awk -v row=$SGE_TASK_ID 'NR == row {gsub(/"/,"",$1); print $1}' "$chrom_set")
af_file="$project_root/$AF_FILE_DIR/TGP_chr${chrom}_afs.txt"

# run the indel script
script="$script_dir/get_common_vars/exon_disruption/get_common_vars.py"
python3 "$script" --af_limit "$AF_LIMIT" --chrom "$chrom" --af_file "$af_file" --exon_file "$project_root/$exon_file" --output_dir "$output_dir" --total_num_chroms "$total_num_chroms"
