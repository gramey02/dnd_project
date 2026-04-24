#!/bin/bash
#$ -N het_combos
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
genes_w_guides="$3"
filtered_vcf_dir="$4"
exon_file="$5"
valid_pairs_fp="$6"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $genes_w_guides)
script="$script_dir/get_hets/het_combos.py"

python3 "$script" --output_dir "$output_dir" \
    --filtered_vcf_dir "$filtered_vcf_dir" \
    --exon_file "$project_root/$exon_file" \
    --gene "$gene" \
    --num_samples "$NUM_SAMPLES" \
    --excise_entire_gene "$EXCISE_ENTIRE_GENE" \
    --valid_pairs_fp "$valid_pairs_fp"
