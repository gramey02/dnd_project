#!/bin/bash
#$ -N excision_guides
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

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

# merge the information on genes targetable by this strategy into one file
BASE="$output_dir"
merged_fp="$output_dir/summary_files/cross_strat_gRNAs/metadata/merged_genes_w_valid_guides_ex.txt"

> "$merged_fp"  # truncate/create file



for strat in excision; do
  fp="$BASE/$strat/excavate/het_individuals/metadata/genes_w_valid_guides.txt"

  if [[ -f "$fp" ]]; then
      awk -v s="$strat" '{print $0 "\t" s}' "$fp" >> "$merged_fp"
  else
      echo "Skipping $strat: file not found ($fp)"
  fi
done

# get the unique genes from the merged file
unique_genes_file="$output_dir/summary_files/cross_strat_gRNAs/metadata/unique_genes_with_valid_guides_excision.txt"
cut -f1 "$merged_fp" | sort -u > "$unique_genes_file"
num_unique_genes=$(wc -l < "$unique_genes_file")
# valid excision pairs
valid_pairs_fp="$output_dir/excision/CommonVars/valid_snp_pairs"
# set script name
shell_script="$script_dir/get_guides/excision_guides_array_setup.sh"

if [[ "$num_unique_genes" -ne 0 ]]; then
  qsub -sync y -t 1-"$num_unique_genes" \
    -l mem_free=5G -l h_rt=24:00:00 \
    -o "$project_root/logs/out/excision_array.out" \
    -e "$project_root/logs/err/excision_array.err" \
    "$shell_script" "$output_dir" "$unique_genes_file" "$param_file" "$valid_pairs_fp"

  echo "Finished excision guide analysis."
else
  echo "No genes with valid guides for excision found"
fi