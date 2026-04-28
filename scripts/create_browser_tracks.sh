#!/bin/bash
#$ -N create_browser_tracks
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o logs/out/create_browser_tracks.out
#$ -e logs/err/create_browser_tracks.err

# load initial variables
project_root="/wynton/home/capra/gramey02/dnd_project"
script_dir="$project_root/scripts"
param_file="$project_root/data/params/params.txt"
source "$param_file"
output_dir="${OUTPUT_DIR}${RUN_NAME}"

# create necessary filepaths
mkdir_if_missing() {
  local d="$1"
  mkdir -p "$d"
}
mkdir_if_missing "$output_dir/summary_files/browser_tracks"
mkdir_if_missing "$output_dir/summary_files/browser_tracks/metadata"
mkdir_if_missing "$output_dir/summary_files/browser_tracks/per_gene_files"
mkdir_if_missing "$output_dir/summary_files/browser_tracks/all_genes_together"

# -------------------------------
# Run BED creation
# -------------------------------
# set script name
bed_script="$script_dir/browser_tracks/create_gene_beds.sh"
# get number of targetable genes across all strategies
merged_fp="$output_dir/summary_files/browser_tracks/metadata/merged_targetable_genes_prepam.txt" # create filename
BASE="$output_dir"
# clear existing file first
merged_fp="$output_dir/summary_files/browser_tracks/metadata/merged_targetable_genes_prepam.txt"
: > "$merged_fp"
for strat in "${EDIT_STRATS[@]}"; do
     # we want to use all genes with common variants, not necessarily NGG-filtered variants
    fp="$BASE/$strat/excavate/input_metadata/excavate_run_metadata.txt"
    [[ -f "$fp" ]] || continue # skips if the file doesn't exist
    awk -v s="$strat" '{print $0 "\t" s}' "$fp" >> "$merged_fp"
done
# get the unique genes from the file
unique_genes_file="$output_dir/summary_files/browser_tracks/metadata/targetable_genes.txt"
cut -f1 "$merged_fp" | sort -u > "$unique_genes_file"
num_targetable_genes=$(wc -l < "$unique_genes_file")
# run the bed creation script
if (( num_targetable_genes == 0 )); then
    echo "No targetable genes found in $unique_genes_file" >&2
else
     qsub -cwd -sync y \
     -l mem_free=5G -l h_rt=02:00:00 \
     -t 1-"$num_targetable_genes" \
     -o "$project_root/logs/out/bed_creation.out" \
     -e "$project_root/logs/err/bed_creation.err" \
     "$bed_script" "$param_file" "$unique_genes_file" "$output_dir"
fi

echo "Finished promoter BED creation."


# --------------------------------------------------------
# Merge resulting files into one bed file across all genes
# --------------------------------------------------------
bed_files=( "$output_dir"/summary_files/browser_tracks/per_gene_files/*/*.bed )

if [[ -e "${bed_files[0]}" ]]; then
    cat "${bed_files[@]}" > "$output_dir/summary_files/browser_tracks/all_genes_together/all_genes.bed"
else
    echo "No BED files found under per_gene_files"
fi


# # -------------------------------
# # Run bigBed conversion
# # -------------------------------
# bigbed_script="$script_dir/browser_tracks/create_bigbeds.sh"

# qsub -cwd \
#      -l mem_free=5G -l h_rt=02:00:00 \
#      -o "$project_root/logs/out/bigbed_conversion.out" \
#      -e "$project_root/logs/err/bigbed_conversion.err" \
#      "$bigbed_script"

# echo "Finished bigBed conversion."