#!/bin/bash
#$ -N excision_guides_dcc
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/exicision_guides_dcc.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/excision_guides_dcc.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

# all genes uncomment
# # input vars
# output_dir=$1
# unique_genes_file=$2
# param_file=$3
# source $param_file
# valid_pairs_fp=$4

# # set up per job gene
# gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $unique_genes_file)
# vcf_dir="$OUTPUT_DIR$RUN_NAME/excision/excavate/Guide_filtered_vcfs" # "$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/downsampled_vcfs"

# OUTPUT_FILE="$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/results/${gene}_excision_gRNAs.csv"

# if [ -f "$OUTPUT_FILE" ]; then
#     echo "Output file for ${gene} exists."
# else
#     # downsample to the top 50 snp locations with the higest heterozygote frequency
#     echo "Output file for ${gene} does not exist, restarting gRNA priortization from scratch or from checkpoint."
#     gRNA_script="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/gRNA_prioritization/excision_guides.py"
#     python3 $gRNA_script --output_dir $output_dir --gene $gene --num_samples $NUM_SAMPLES --valid_pairs_fp $valid_pairs_fp --vcf_dir "$vcf_dir" --max_iter 50
# fi

# ---------------------------------------------------------------
# single fallout gene below
# input vars
param_file="$project_root/data/params/params.txt"
source "$param_file"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

valid_pairs_fp="$resolved_output_base$RUN_NAME/excision/CommonVars/valid_snp_pairs"

# set up per job gene
gene='DCC'
vcf_dir="$resolved_output_base$RUN_NAME/excision/excavate/Guide_filtered_vcfs" # "$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/downsampled_vcfs"
output_dir="$resolved_output_base$RUN_NAME"

OUTPUT_FILE="$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/results/${gene}_excision_gRNAs.csv"

if [ -f "$OUTPUT_FILE" ]; then
    echo "Output file for ${gene} exists."
else
    # downsample to the top 258 snp locations with the higest heterozygote frequency
    echo "Output file for ${gene} does not exist, restarting gRNA priortization from scratch or from checkpoint."
    gRNA_script="$script_dir/excision_guides.py"
    python3 "$gRNA_script" --output_dir "$output_dir" --gene "$gene" --num_samples "$NUM_SAMPLES" --valid_pairs_fp "$valid_pairs_fp" --vcf_dir "$vcf_dir" --max_iter 50
fi
