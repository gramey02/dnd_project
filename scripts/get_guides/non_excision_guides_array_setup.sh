#!/bin/bash
#$ -N non_excision_guides
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/non_exicision_guides.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/non_excision_guides.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# input vars
output_dir="$1"
unique_genes_file="$2"
param_file="$3"
source "$param_file"
all_strats_together="$4"

# set up per job gene
gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' "$unique_genes_file")
echo "$gene"

gRNA_script="$script_dir/non_excision_guides.py"
python3 "$gRNA_script" --output_dir "$output_dir" --gene "$gene" --num_samples "$NUM_SAMPLES" --all_strats_together "$all_strats_together"
