#!/bin/bash
#$ -N non_excision_guides
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ../../logs/out/non_excision_guides.out
#$ -e ../../logs/err/non_excision_guides.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"
execution_utils="$project_root/scripts/utils/execution_mode.sh"

output_dir="$1"
param_file="$2"
source "$param_file"
source "$execution_utils"

# make directory to hold info
mkdir_if_missing() {
  local d="$1"
  mkdir -p "$d"
}
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/metadata"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/checkpoints"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/results"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/logs"

# merge the information on genes into one file
BASE="$output_dir"
merged_fp="$output_dir/summary_files/cross_strat_gRNAs/metadata/merged_genes_w_valid_guides.txt"

> "$merged_fp"  # truncate/create file

for strat in indels CRISPRoff acceptor_base_edits donor_base_edits; do
    awk -v s="$strat" '{print $0 "\t" s}' \
        "$BASE/$strat/excavate/het_individuals/metadata/genes_w_valid_guides.txt" \
        >> "$merged_fp"
done

# get the unique genes from the file
unique_genes_file="$output_dir/summary_files/cross_strat_gRNAs/metadata/unique_genes_with_valid_guides_non_excision.txt"
cut -f1 "$merged_fp" | sort -u > "$unique_genes_file"
num_unique_genes=$(wc -l < "$unique_genes_file")

# set up variable to run all strategies together or separately
all_strats_together="False"

shell_script="$script_dir/non_excision_guides_array_setup.sh"
qsub -t 1-"$num_unique_genes" -l mem_free=5G -l h_rt=5:00:00 $shell_script "$output_dir" "$unique_genes_file" "$param_file" $all_strats_together
