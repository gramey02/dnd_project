#!/bin/bash
#$ -N non_excision_guides
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input variables
output_dir="$1"
param_file="$2"
source "$param_file"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# echo "Project root"
# echo $project_root
# echo "Script dir"
# echo $script_dir

# make directory to hold guide outputs
mkdir_if_missing() {
  local d="$1"
  mkdir -p "$d"
}
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/non_excision_guides"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/non_excision_guides/metadata"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/non_excision_guides/checkpoints"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/non_excision_guides/results"
mkdir_if_missing "$output_dir/summary_files/cross_strat_gRNAs/non_excision_guides/logs"

# merge the information on genes targetable by these strategies into one file
BASE="$output_dir"
merged_fp="$output_dir/summary_files/cross_strat_gRNAs/non_excision_guides/metadata/merged_genes_w_valid_guides_nonex.txt"

> "$merged_fp"  # truncate/create file

# echo "Base"
# echo $BASE
# echo "Merged filepath"
# echo $merged_fp

# create a list of the genes with guides across all non-excision strategies
for strat in indels CRISPRoff acceptor_base_edits donor_base_edits; do
    fp="$BASE/$strat/excavate/het_individuals/metadata/genes_w_valid_guides.txt"

    [[ -f "$fp" ]] || continue # skips if the file doesn't exist

    awk -v s="$strat" '{print $0 "\t" s}' "$fp" >> "$merged_fp"
done

# get the unique genes from the file
unique_genes_file="$output_dir/summary_files/cross_strat_gRNAs/non_excision_guides/metadata/unique_genes_with_valid_guides_non_excision.txt"
cut -f1 "$merged_fp" | sort -u > "$unique_genes_file"
num_unique_genes=$(wc -l < "$unique_genes_file")

# echo "Unique genes file"
# echo $unique_genes_file

# set up variable to run all strategies together or separately
all_strats_together="False"
# set shell script name
shell_script="$script_dir/get_guides/non_excision_guides_array_setup.sh"

# run
if [[ "$num_unique_genes" -ne 0 ]]; then
  qsub -sync y -t 1-"$num_unique_genes" \
    -l mem_free=5G -l h_rt=05:00:00 \
    -o "$project_root/logs/out/non_excision_array.out" \
    -e "$project_root/logs/err/non_excision_array.err" \
    "$shell_script" "$output_dir" "$unique_genes_file" "$param_file" "$all_strats_together"

  echo "Finished non-excision guide analysis."
else
  echo "No genes with valid guides for non-excision found"
fi