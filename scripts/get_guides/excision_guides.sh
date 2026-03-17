#!/bin/bash
#$ -N excision_setup
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/exicision_setup2.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/excision_setup2.err

param_file="/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/params.txt"
source $param_file

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
mkdir_if_missing "$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs"
mkdir_if_missing "$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/metadata"
mkdir_if_missing "$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/checkpoints"
mkdir_if_missing "$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/results"
mkdir_if_missing "$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/logs"

# merge the information on genes into one file
BASE=$OUTPUT_DIR$RUN_NAME
merged_fp="$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/metadata/merged_genes_w_valid_guides.txt"

> "$merged_fp"  # truncate/create file

for strat in excision; do
    awk -v s="$strat" '{print $0 "\t" s}' \
        "$BASE/$strat/excavate/het_individuals/metadata/genes_w_valid_guides.txt" \
        >> "$merged_fp"
done

# get the unique genes from the file
unique_genes_file="$OUTPUT_DIR$RUN_NAME/summary_files/cross_strat_gRNAs/metadata/unique_genes_with_valid_guides_non_excision.txt"
cut -f1 "$merged_fp" | sort -u > "$unique_genes_file"
num_unique_genes=$(wc -l < "$unique_genes_file")
# valid excision pairs
valid_pairs_fp="$OUTPUT_DIR$RUN_NAME/excision/CommonVars/valid_snp_pairs"

# run array job prioritizing non-excision gRNAs for each gene
shell_script="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/gRNA_prioritization/excision_guides_array_setup.sh"
#gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $unique_genes_file)
output_dir="$OUTPUT_DIR$RUN_NAME"
qsub -t 1-"$num_unique_genes" -l mem_free=5G -l h_rt=24:00:00 $shell_script "$output_dir" "$unique_genes_file" "$param_file" "$valid_pairs_fp"
