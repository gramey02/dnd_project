#!/bin/bash
#$ -N het_combos
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ../../logs/out/het_combos.out
#$ -e ../../logs/err/het_combos.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
genes_w_guides="$3"
filtered_vcf_dir="$4"
exon_file="$5"

gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $genes_w_guides)
script="$script_dir/het_combos.py"

python3 "$script" --output_dir "$output_dir" \
    --filtered_vcf_dir "$filtered_vcf_dir" \
    --exon_file "$exon_file" \
    --gene_info "$genes_w_guides" \
    --gene "$gene" \
    --num_samples "$NUM_SAMPLES" \
    --excise_entire_gene "$EXCISE_ENTIRE_GENE"
