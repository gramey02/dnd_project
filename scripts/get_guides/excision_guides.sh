#!/bin/bash
#$ -N excision_setup
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/exicision_setup2.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/excision_setup2.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

param_file="$project_root/data/params/params.txt"
source "$param_file"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

# make directory to hold info
mkdir_if_missing() {
  local d="$1"
  if [ -d "$d" ]; then
    echo "Directory '$d' already exists."
  else
    mkdir -p "$d"       # -p creates parents as needed; safe to call repeatedly
    echo "Directory '$d' created."
  fi
}
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs"
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/metadata"
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/checkpoints"
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/results"
mkdir_if_missing "$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/logs"

# merge the information on genes into one file
BASE="$resolved_output_base$RUN_NAME"
merged_fp="$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/metadata/merged_genes_w_valid_guides.txt"

> "$merged_fp"  # truncate/create file

for strat in excision; do
    awk -v s="$strat" '{print $0 "\t" s}' \
        "$BASE/$strat/excavate/het_individuals/metadata/genes_w_valid_guides.txt" \
        >> "$merged_fp"
done

# get the unique genes from the file
unique_genes_file="$resolved_output_base$RUN_NAME/summary_files/cross_strat_gRNAs/metadata/unique_genes_with_valid_guides_non_excision.txt"
cut -f1 "$merged_fp" | sort -u > "$unique_genes_file"
num_unique_genes=$(wc -l < "$unique_genes_file")
# valid excision pairs
valid_pairs_fp="$resolved_output_base$RUN_NAME/excision/CommonVars/valid_snp_pairs"

# run array job prioritizing non-excision gRNAs for each gene
shell_script="$script_dir/excision_guides_array_setup.sh"
#gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $unique_genes_file)
output_dir="$resolved_output_base$RUN_NAME"
qsub -t 1-"$num_unique_genes" -l mem_free=5G -l h_rt=24:00:00 "$shell_script" "$output_dir" "$unique_genes_file" "$param_file" "$valid_pairs_fp"
