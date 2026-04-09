#!/bin/bash
#$ -N excision_guides_array_setup
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input args
output_dir="$1"
unique_genes_file="$2"
param_file="$3"
source "$param_file"
valid_pairs_fp="$4"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"


# Resolve the one gene assigned to this local loop iteration.
gene=$(awk -v row="$SGE_TASK_ID" 'NR == row {print $1}' "$unique_genes_file")
vcf_dir="$OUTPUT_DIR$RUN_NAME/excision/excavate/Guide_filtered_vcfs" # "$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/downsampled_vcfs"
OUTPUT_FILE="$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/results/${gene}_excision_gRNAs.csv"

# run greedy algorithm on each gene, with max 50 iterations
if [ ! -f "$OUTPUT_FILE" ]; then
    # downsample to the top 258 snp locations with the higest heterozygote frequency
    gRNA_script="$script_dir/scripts/get_guides/excision_guides.py"
    python3 "$gRNA_script" --output_dir "$output_dir" --gene "$gene" --num_samples "$NUM_SAMPLES" --valid_pairs_fp "$valid_pairs_fp" --vcf_dir "$vcf_dir" --max_iter 50
fi
