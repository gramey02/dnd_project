#!/bin/bash
#$ -N crisproff_pipeline
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input arguments
param_file="$2"
source "$param_file"
output_dir="$1"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# find exon_file
exon_file="$EXON_FILE_FOR_ANALYSIS"

# get num chromosomes
if [ ! -e "filename.txt" ]; then
    col=$(head -1 "$project_root/$exon_file" \
    | tr ',' '\n' \
    | sed 's/"//g' \
    | grep -nx "chromosome_name" \
    | cut -d: -f1)
    # fail early if column not found
    if [[ -z "$col" ]]; then
        echo "ERROR: chromosome_name column not found in $project_root/$exon_file"
        exit 1
    fi

    # save chromosomes into their own file
    cut -d',' -f"$col" "$project_root/$exon_file" | tail -n +2 | sort -u > "$OUTPUT_DIR/$RUN_NAME/chromosomes/chrom_set.txt"
    num_chroms=$(wc -l < "$OUTPUT_DIR/$RUN_NAME/chromosomes/chrom_set.txt")
fi


# script to get ubiquitous promoter regions & identify common vars in them
promoter_common_vars="$project_root/scripts/get_common_vars/epi_silencing/promoter_common_vars.sh"
echo "Started identifying ubiquitous promoter regions & common vars..."
qsub -sync y -l mem_free=1G -l h_rt=00:20:00 -o "$project_root/logs/out/promoter_common_vars.out" -e "$project_root/logs/err/promoter_common_vars.err" "$promoter_common_vars" "$output_dir" "$param_file" "$exon_file"
echo "Finished identifiying ubiquitous promoter regions & common vars."

# get number of genes that have common vars in promoter regions
common_var_genes="$output_dir/ubiq_region_CommonVars/CommonVars_VarNumOver0_summary_noIDX.txt"
awk -F'\t' '$2+0 > 0' "$output_dir/ubiq_region_CommonVars/CommonVars_ALL_summary_noIDX.txt" > "$common_var_genes"
num_common_var_genes=$(wc -l < "$common_var_genes") # get the number of genes that have common vars in them (from summary file)

# script to generate text files that we'll use to filter vcfs (which we'll then input into EXCAVATE)
generate_variant_textFiles="$project_root/scripts/format_variants/generate_variant_textFiles.py"
echo "Started generating common var loc files..."
cv_dict_filepath="$output_dir/ubiq_region_CommonVars/CommonVars_ALL_dict.pkl"
python3 "$generate_variant_textFiles" --cv_dict_filepath "$cv_dict_filepath" --exon_file "$project_root/$exon_file" --output_dir "$output_dir" --af_file_dir "$project_root/$AF_FILE_DIR"
echo "Finished generating common var loc files."

# script to filter vcfs accordingly
excavate_vcf_creation="$project_root/scripts/format_variants/generate_filtered_vcfs.sh" # first we need to generate vcf.gz files for the genes that have variants in their ubiquitous regions
input_metadata="$output_dir/excavate/input_metadata/excavate_run_metadata.txt"
echo "Started creating vcf files for excavate input..."
qsub -t 1-"$num_common_var_genes" -l mem_free=2G -l h_rt=01:00:00 -sync y -o "$project_root/logs/out/filt_vcfs_crisproff.out" -e "$project_root/logs/err/filt_vcfs_crisproff.err" "$excavate_vcf_creation" "$output_dir" "$param_file" "$input_metadata"
echo "Finished creating excavate inputs."

# script to run excavate
run_excavate_script="$project_root/scripts/excavate/run_excavate.sh"
echo "Started running EXCAVATE..."
qsub -t 1-"$num_common_var_genes" -l mem_free=1G -l h_rt=00:45:00 -sync y -o "$project_root/logs/out/excavate_crisproff.out" -e "$project_root/logs/err/excavate_crisproff.err" "$run_excavate_script" "$output_dir" "$param_file" "$input_metadata"
echo "Finished running EXCAVATE."

# generate text files for the valid guides so you can filter the vcfs accordingly
generate_guide_textFiles="$project_root/scripts/format_variants/generate_guide_textFiles.py"
guides_filepath="$output_dir/excavate/excavate_outputs" # switch to not include indels later
echo "Started generating position files for valid guides..."
python3 "$generate_guide_textFiles" --guides_filepath "$guides_filepath" --exon_file "$project_root/$exon_file" --output_dir "$output_dir"
echo "Finished generating position files for valid guides."

# filter vcfs based on these viable EXCAVATE guides to see how many heterozygous individuals they will capture
guide_based_filtering="$project_root/scripts/format_variants/position_filtering.sh"
genes_w_guides="$output_dir/excavate/het_individuals/metadata/genes_w_valid_guides.txt"
num_genes_w_guides=$(awk -F'\t' '$1 != "" {n++} END{print n}' "$genes_w_guides")
echo "Started filtering vcfs based on valid guides..."
qsub -t 1-"$num_genes_w_guides" -l mem_free=2G -l h_rt=01:00:00 -sync y -o "$project_root/logs/out/guide_filtering_crisproff.out" -e "$project_root/logs/err/guide_filtering_crisproff.err" "$guide_based_filtering" "$output_dir" "$param_file" "$genes_w_guides"
echo "Finished filtering vcfs based on valid guides."

# script to calculate number of heterozygous individuals for each gene
get_targeted_hets="$project_root/scripts/get_hets/get_targeted_hets.sh"
echo "Started calculating heterozygous indvidual numbers..."
excavate_output_dir="$output_dir/excavate/excavate_outputs"
filtered_vcf_dir="$output_dir/excavate/Guide_filtered_vcfs"
qsub -l mem_free=2G -l h_rt=01:00:00 -sync y -o "$project_root/logs/out/hets_crisproff.out" -e "$project_root/logs/err/hets_crisproff.err" "$get_targeted_hets" "$output_dir/excavate/het_individuals" "$param_file" "$genes_w_guides" "$excavate_output_dir" "$filtered_vcf_dir"
echo "Finished calculating heterozygous individual numbers."
