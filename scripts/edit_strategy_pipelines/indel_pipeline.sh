#!/bin/bash
#$ -N indel_pipeline
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input arguments
param_file="$2"
source "$param_file"
output_dir="$1"
# find exon_file
exon_file="$EXON_FILE_FOR_ANALYSIS"
# set top level directory information
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# # get number of chromosomes represented in exon file
# #col=$(head -1 "$project_root/$exon_file" | tr ',' '\n' | grep -nx "chromosome_name" | cut -d: -f1)
# col=$(head -1 "$project_root/$exon_file" \
#     | tr ',' '\n' \
#     | sed 's/"//g' \
#     | grep -nx "chromosome_name" \
#     | cut -d: -f1)
# # fail early if column not found
# if [[ -z "$col" ]]; then
#     echo "ERROR: chromosome_name column not found in $project_root/$exon_file"
#     exit 1
# fi
# # save chromosomes into their own file
# cut -d',' -f"$col" "$project_root/$exon_file" | tail -n +2 | sort -u > "$OUTPUT_DIR/$RUN_NAME/chromosomes/chrom_set.txt"
# num_chroms=$(wc -l < "$OUTPUT_DIR/$RUN_NAME/chromosomes/chrom_set.txt")

# # script to get ubiquitous exonic regions & identify common vars in them
# ubiq_regions_CommonVar_script="$project_root/scripts/get_common_vars/exon_disruption/get_common_vars_UER.sh"
# echo "Started identifying common coding vars..."
# qsub -l mem_free=1G -l h_rt=04:20:00 -t 1-"$num_chroms" -sync y -o "$project_root/logs/out/get_common_vars_UER.out" -e "$project_root/logs/err/get_common_vars_UER.err" "$ubiq_regions_CommonVar_script" "$output_dir" "$param_file" "$num_chroms" "$exon_file"
# echo "Finished identifiying common coding vars."

# # merge per-chromosome common var outputs (single serial step, run once the array job has finished)
# merge_common_vars_UER="$project_root/scripts/get_common_vars/exon_disruption/merge_common_vars_UER.sh"
# echo "Merging per-chromosome common coding var outputs..."
# qsub -sync y -l mem_free=2G -l h_rt=00:15:00 -o "$project_root/logs/out/merge_common_vars_UER.out" -e "$project_root/logs/err/merge_common_vars_UER.err" "$merge_common_vars_UER" "$output_dir" "$param_file" "$num_chroms"
# if [ $? -ne 0 ]; then
#     echo "ERROR: merge_common_vars_UER job failed (see logs/err/merge_common_vars_UER.err)"
#     exit 1
# fi
# echo "Finished merging per-chromosome common coding var outputs."

# ##### NMD script goes here (only runs when RUN_NMD is enabled)
# targetable_common_var_dict="$output_dir/ubiq_region_CommonVars/CommonVars_ALL_dict.pkl"
# if [[ "$RUN_NMD" == "1" || "$RUN_NMD" == "True" || "$RUN_NMD" == "true" ]]; then
#     NMD_script="$project_root/scripts/NMD/Annot_NMD_escape.sh"
#     echo "Started running NMD analysis on guide-filtered indel vars..."
#     qsub -l mem_free=1G -l h_rt=05:00:00 -sync y -o "$project_root/logs/out/NMD.out" -e "$project_root/logs/err/NMD.err" "$NMD_script" "$output_dir" "$param_file" "$targetable_common_var_dict"
#     echo "Finished running NMD analysis."
# else
#     # NMD disabled: format common vars the same way the other edit-strategy pipelines do
#     generate_variant_textFiles="$project_root/scripts/format_variants/generate_variant_textFiles.py"
#     echo "Started generating common var loc files..."
#     python3 "$generate_variant_textFiles" --cv_dict_filepath "$targetable_common_var_dict" --exon_file "$project_root/$exon_file" --output_dir "$output_dir" --af_file_dir "$project_root/$AF_FILE_DIR"
#     echo "Finished generating common var loc files."
# fi
common_var_genes="$output_dir/excavate/input_metadata/excavate_run_metadata.txt"
num_common_var_genes=$(wc -l < "$common_var_genes") # get the number of genes that have common vars in them (from summary file)

# if [[ "$num_common_var_genes" -eq 0 ]]; then
#     echo "No genes with common vars in shared exons found. Nothing to run for indel strategy — exiting."
#     exit 0
# fi

# # script to filter vcfs accordingly
# excavate_vcf_creation="$project_root/scripts/format_variants/generate_filtered_vcfs.sh" # first we need to generate vcf.gz files for the genes that have variants in their ubiquitous exonic regions
input_metadata="$output_dir/excavate/input_metadata/excavate_run_metadata.txt"
# echo "Started creating vcf files for excavate input..."
# qsub -t 1-"$num_common_var_genes" -l mem_free=2G -l h_rt=05:00:00 -sync y -o "$project_root/logs/out/filt_vcfs_indels.out" -e "$project_root/logs/err/filt_vcfs_indels.err" "$excavate_vcf_creation" "$output_dir" "$param_file" "$input_metadata"
# echo "Finished creating excavate inputs."

# script to run excavate
# run_excavate_script="$project_root/scripts/excavate/run_excavate.sh"
# echo "Started running EXCAVATE..."
# qsub -t 1-"$num_common_var_genes" -l mem_free=4G -l h_rt=04:45:00 -sync y -o "$project_root/logs/out/excavate_indels.out" -e "$project_root/logs/err/excavate_indels.err" "$run_excavate_script" "$output_dir" "$param_file" "$input_metadata"
# echo "Finished running EXCAVATE."

# one-off run: excavate with off-targets, saved to excavate/excavate_outputs_w_offtargets
run_excavate_script="$project_root/scripts/excavate/run_excavate_w_offtargets.sh"
echo "Started running EXCAVATE (with off-targets)..."
qsub -t 1-"$num_common_var_genes" -l mem_free=20G -l h_rt=08:00:00 -sync y -o "$project_root/logs/out/excavate_indels_w_offtargets.out" -e "$project_root/logs/err/excavate_indels_w_offtargets.err" "$run_excavate_script" "$output_dir" "$param_file" "$input_metadata"
echo "Finished running EXCAVATE (with off-targets)."

# # generate text files for the valid guides so you can filter the vcfs accordingly
# generate_guide_textFiles="$project_root/scripts/format_variants/generate_guide_textFiles.py"
# guides_filepath="$output_dir/excavate/excavate_outputs" # switch to not include indels later
# echo "Started generating position files for valid guides..."
# python3 "$generate_guide_textFiles" --guides_filepath "$guides_filepath" --exon_file "$project_root/$exon_file" --output_dir "$output_dir"
# echo "Finished generating position files for valid guides."

# # filter vcfs based on these viable EXCAVATE guides to see how many heterozygous individuals they will capture
# guide_based_filtering="$project_root/scripts/format_variants/position_filtering.sh"
# genes_w_guides="$output_dir/excavate/het_individuals/metadata/genes_w_valid_guides.txt"
# num_genes_w_guides=$(awk -F'\t' '$1 != "" {n++} END{print n}' "$genes_w_guides")
# if [[ "$num_genes_w_guides" -gt 0 ]]; then
#     echo "Started filtering vcfs based on valid guides..."
#     qsub -t 1-"$num_genes_w_guides" -l mem_free=2G -l h_rt=05:00:00 -sync y -o "$project_root/logs/out/guide_filtering_indels.out" -e "$project_root/logs/err/guide_filtering_indels.err" "$guide_based_filtering" "$output_dir" "$param_file" "$genes_w_guides"
#     echo "Finished filtering vcfs based on valid guides."
# else
#     echo "No genes with valid guides found. Skipping position filtering."
# fi

# # script to calculate the number of heterozygous individuals for each gene, pre-excavate filtering
# get_targeted_hets_prePAM="$project_root/scripts/get_hets/get_targeted_hets_prePAM.sh"
# echo "Started calculating heterozygous individual numbers (prePAM filtering)..."
# filtered_vcf_dir="$output_dir/excavate/input_vcfs"
# gene_info="$output_dir/excavate/input_metadata/excavate_run_metadata.txt"
# qsub -l mem_free=2G -l h_rt=04:00:00 -o "$project_root/logs/out/hets_indels_prepam.out" -e "$project_root/logs/err/hets_indels_prepam.err" "$get_targeted_hets_prePAM" "$output_dir/prePAM_hets" "$param_file" "$gene_info" "$filtered_vcf_dir"

# # script to calculate number of heterozygous individuals for each gene
# if [[ "$num_genes_w_guides" -gt 0 ]]; then
#     get_targeted_hets="$project_root/scripts/get_hets/get_targeted_hets.sh"
#     echo "Started calculating heterozygous indvidual numbers..."
#     excavate_output_dir="$output_dir/excavate/excavate_outputs"
#     filtered_vcf_dir="$output_dir/excavate/Guide_filtered_vcfs"
#     qsub -l mem_free=2G -l h_rt=04:00:00 -sync y -o "$project_root/logs/out/hets_indels.out" -e "$project_root/logs/err/hets_indels.err" "$get_targeted_hets" "$output_dir/excavate/het_individuals" "$param_file" "$genes_w_guides" "$excavate_output_dir" "$filtered_vcf_dir"
#     echo "Finished calculating heterozygous individual numbers."
# else
#     echo "No genes with valid guides found. Skipping heterozygous-individual calculation."
# fi
