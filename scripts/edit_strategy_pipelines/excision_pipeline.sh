#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

# parse input arguments
param_file="$2"
source "$param_file"
output_dir="$1"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

# Local replacement for SGE array jobs: call the target script once for each
# 1-based row index that would previously have come from SGE_TASK_ID.
run_indexed_jobs() {
  local count="$1"
  shift
  if (( count <= 0 )); then
    return 0
  fi

  local task_id
  for ((task_id=1; task_id<=count; task_id++)); do
    bash "$@" "$task_id"
  done
}

# set exon file
#exon_file=$EXON_FILE_FOR_ANALYSIS
# exon_file="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/filtered_exon_info/dHS_exonInfo_GTExCurated_transcripts_TPM1.0_ExprProp0.01_2025_04_22.csv"
#exon_file="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/filtered_exon_info/donor_test_exon.csv"
exon_file="$EXON_FILE_FOR_ANALYSIS"

# script to generate appropriate coordinates to look for excision variants in
find_variant_region="$project_root/scripts/format_variants/find_variant_region.py"
python3 "$find_variant_region" --exon_file "$exon_file" --output_dir "$output_dir/CommonVars" --upstream_excision_dist "$UPSTREAM_EXCISION_DIST" --downstream_excision_dist "$DOWNSTREAM_EXCISION_DIST" --nearest_gene_file "$GENE_BODY_FILE" --filter_out_nearby_genes "$FILTER_OUT_NEARBY_GENES"

# script to generate list of variants in excision windows to run through excavate
find_commonVars="$project_root/scripts/get_common_vars/excision/find_excision_commonVars.sh"
gene_info="$output_dir/CommonVars/gene_filtered_excision_coords.txt"
echo "Started extracting common var info..."
bash "$find_commonVars" "$output_dir/CommonVars" "$param_file" "$gene_info" "$exon_file"
echo "Finished extracting common variant info."

# get number of genes that have common vars in excision_window --> make sure there are 2 or more, since that's the minimum required for excision
common_var_genes=$output_dir"/CommonVars/CommonVars_VarNumOver2_summary_noIDX.txt"
awk -F'\t' '$2+0 >= 2' "$output_dir/CommonVars/CommonVars_ALL_summary_noIDX.txt" > "$common_var_genes"
num_common_var_genes=$(wc -l < $common_var_genes) # get the number of genes that have common vars in them (from summary file)
# I think the filtering in the next step is better

# script to determine viable pairs of snps, to narrow down the number of snps for future steps
filter_excision_snps="$project_root/scripts/format_variants/filter_excision_snps.sh"
cv_dict_filepath=$output_dir"/CommonVars/CommonVars_ALL_dict.pkl"
echo "Filtering snps based on those that encompass exons across transcripts..."
run_indexed_jobs "$num_common_var_genes" "$filter_excision_snps" "$output_dir/CommonVars/refined_common_vars" "$param_file" "$cv_dict_filepath" "$exon_file" "$common_var_genes"
echo "Finished filtering excision snps."

# script to generate text files that we'll use to filter vcfs (which we'll then input into EXCAVATE)
generate_variant_textFiles="$project_root/scripts/format_variants/generate_variant_textFiles.py"
echo "Started generating common var loc files..."
cv_dict_filepath=$output_dir"/CommonVars/CommonVars_ALL_dict.pkl"
python3 "$generate_variant_textFiles" --cv_dict_filepath "$cv_dict_filepath" --exon_file "$exon_file" --output_dir "$output_dir" --af_file_dir "$AF_FILE_DIR" --edit_strat "excision"
num_common_var_genes=$(ls -1 $output_dir"/excavate/CommonVar_locs" | wc -l)
echo "Finished generating common var loc files."

# script to filter vcfs accordingly
excavate_vcf_creation="$project_root/scripts/format_variants/generate_filtered_vcfs.sh" # first we need to generate vcf.gz files for the genes that have variants in their ubiquitous regions
input_metadata=$output_dir"/excavate/input_metadata/excavate_run_metadata.txt"
echo "Started creating vcf files for excavate input..."
run_indexed_jobs "$num_common_var_genes" "$excavate_vcf_creation" "$output_dir" "$param_file" "$input_metadata"
echo "Finished creating excavate inputs."

# script to run excavate
run_excavate_script="$project_root/scripts/EXCAVATE_HT_run/run_excavate.sh"
echo "Started running EXCAVATE..."
run_indexed_jobs "$num_common_var_genes" "$run_excavate_script" "$output_dir" "$param_file" "$input_metadata"
echo "Finished running EXCAVATE."

# generate text files for the valid guides so you can filter the individuals' vcfs accordingly
generate_guide_textFiles="$project_root/scripts/format_variants/generate_guide_textFiles.py"
guides_filepath=$output_dir"/excavate/excavate_outputs"
echo "Started generating position files for valid guides..."
python3 "$generate_guide_textFiles" --guides_filepath "$guides_filepath" --exon_file "$exon_file" --output_dir "$output_dir"
echo "Finished generating position files for valid guides."

# filter vcfs based on these viable EXCAVATE guides to see how many heterozygous individuals they will capture
guide_based_filtering="$project_root/scripts/format_variants/position_filtering.sh"
genes_w_guides=$output_dir"/excavate/het_individuals/metadata/genes_w_valid_guides.txt"
num_genes_w_guides=$(awk -F'\t' '$1 != "" {n++} END{print n}' $genes_w_guides)
echo "Started filtering vcfs based on valid guides..."
run_indexed_jobs "$num_genes_w_guides" "$guide_based_filtering" "$output_dir" "$param_file" "$genes_w_guides"
echo "Finished filtering vcfs based on valid guides."

# script to calculate number of heterozygous individuals that will be hit by pairs of snps (also runs greedy algorithm)
het_combos_script="$project_root/scripts/get_hets/het_combos.sh"
filtered_vcf_dir=$output_dir"/excavate/Guide_filtered_vcfs"
echo "Started capturing heterozygote excision information..."
echo $exon_file
run_indexed_jobs "$num_genes_w_guides" "$het_combos_script" "$output_dir/excavate/het_individuals" "$param_file" "$genes_w_guides" "$filtered_vcf_dir" "$exon_file"
echo "Finished capturing heterozygote excision information."

# script to get the number of individuals haplotype-specific guides per gene
get_guide_info="$project_root/scripts/get_guides/excision_guides.sh"
echo "Calculating number of guides to target heterozygotes..."
bash "$get_guide_info" "$output_dir/excavate/guide_numbers" "$param_file" "$genes_w_guides" "$filtered_vcf_dir" "$exon_file"
echo "Finished calculating number of guides."
