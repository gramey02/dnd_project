#!/bin/bash
#$ -N indel_pipeline
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/indel_pipeline.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/indel_pipeline.err

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

# find exon_file
#exon_file=$EXON_FILE_FOR_ANALYSIS
#exon_file="/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/filtered_exon_info/test_exon.csv"
#exon_file="/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/filtered_exon_info/dHS_exonInfo_GTExCurated_transcripts_TPM1.0_ExprProp0.01_2025_04_22.csv"
exon_file="$EXON_FILE_FOR_ANALYSIS"

# get num chromosomes
col=$(head -1 "$exon_file" | tr ',' '\n' | grep -nx "chromosome_name" | cut -d: -f1)
cut -d',' -f"$col" "$exon_file" | tail -n +2 | sort -u > "$resolved_output_base$RUN_NAME/chromosomes/chrom_set.txt"
num_chroms=$(wc -l < "$resolved_output_base$RUN_NAME/chromosomes/chrom_set.txt")

# script to get ubiquitous exonic regions & identify common vars in them
ubiq_regions_CommonVar_script="$project_root/scripts/get_common_vars/exon_disruption/get_common_vars_UER.sh"
echo "Started identifying common coding vars..."
qsub -l mem_free=1G -l h_rt=00:20:00 -t 1-"$num_chroms" -sync y "$ubiq_regions_CommonVar_script" "$output_dir" "$param_file" "$num_chroms" "$exon_file"
echo "Finished identifiying common coding vars."

##### NMD script goes here
targetable_common_var_dict="$output_dir/ubiq_region_CommonVars/CommonVars_ALL_dict.pkl"
NMD_script="$project_root/scripts/NMD/Annot_NMD_escape.sh"
echo "Started running NMD analysis on guide-filtered indel vars..."
qsub -l mem_free=1G -l h_rt=01:00:00 -sync y "$NMD_script" "$output_dir" "$param_file" "$targetable_common_var_dict"
echo "Finished running NMD analysis."
common_var_genes=$output_dir"/excavate/input_metadata/excavate_run_metadata.txt"
num_common_var_genes=$(wc -l < $common_var_genes) # get the number of genes that have common vars in them (from summary file)

#--------------these are not necessary for indel pipeline, since the NMD pipeline takes care of vcf formatting--------------
# # scripts to generate excavate input files
# common_var_genes=$output_dir"/ubiq_region_CommonVars/CommonVars_VarNumOver0_summary_noIDX.txt"
# awk -F'\t' '$2+0 > 0' $output_dir"/ubiq_region_CommonVars/CommonVars_ALL_summary_noIDX.txt" > $common_var_genes
# num_common_var_genes=$(wc -l < $common_var_genes) # get the number of genes that have common vars in them (from summary file)

# # script to generate text files that we'll use to filter vcfs (which we'll then input into EXCAVATE)
# generate_variant_textFiles="/wynton/protected/home/capra/gramey02/ConklinCollab/scripts/DN_PAMsites/generate_variant_textFiles.py"
# echo "Started generating common var loc files..."
# cv_dict_filepath=$output_dir"/ubiq_region_CommonVars/CommonVars_ALL_dict.pkl"
# python3 $generate_variant_textFiles --cv_dict_filepath $cv_dict_filepath --exon_file $exon_file --output_dir $output_dir
# echo "Finished generating common var loc files."
#------------------------------------------------------------------------------------------------------------------------------

# script to filter vcfs accordingly
excavate_vcf_creation="$project_root/scripts/format_variants/generate_filtered_vcfs.sh" # first we need to generate vcf.gz files for the genes that have variants in their ubiquitous exonic regions
input_metadata=$output_dir"/excavate/input_metadata/excavate_run_metadata.txt"
echo "Started creating vcf files for excavate input..."
qsub -t 1-"$num_common_var_genes" -l mem_free=2G -l h_rt=01:00:00 -sync y "$excavate_vcf_creation" "$output_dir" "$param_file" "$input_metadata"
echo "Finished creating excavate inputs."

# script to run excavate
run_excavate_script="$project_root/scripts/EXCAVATE_HT_run/run_excavate.sh"
echo "Started running EXCAVATE..."
qsub -t 1-"$num_common_var_genes" -l mem_free=1G -l h_rt=00:45:00 -sync y "$run_excavate_script" "$output_dir" "$param_file" "$input_metadata"
echo "Finished running EXCAVATE."

# generate text files for the valid guides so you can filter the vcfs accordingly
generate_guide_textFiles="$project_root/scripts/format_variants/generate_guide_textFiles.py"
guides_filepath=$output_dir"/excavate/excavate_outputs" # switch to not include indels later
echo "Started generating position files for valid guides..."
python3 "$generate_guide_textFiles" --guides_filepath "$guides_filepath" --exon_file "$exon_file" --output_dir "$output_dir"
echo "Finished generating position files for valid guides."

# filter vcfs based on these viable EXCAVATE guides to see how many heterozygous individuals they will capture
guide_based_filtering="$project_root/scripts/format_variants/position_filtering.sh"
genes_w_guides=$output_dir"/excavate/het_individuals/metadata/genes_w_valid_guides.txt"
num_genes_w_guides=$(awk -F'\t' '$1 != "" {n++} END{print n}' $genes_w_guides)
echo "Started filtering vcfs based on valid guides..."
qsub -t 1-"$num_genes_w_guides" -l mem_free=2G -l h_rt=01:00:00 -sync y "$guide_based_filtering" "$output_dir" "$param_file" "$genes_w_guides"
echo "Finished filtering vcfs based on valid guides."

# script to calculate number of heterozygous individuals for each gene
get_targeted_hets="$project_root/scripts/get_hets/get_targeted_hets.sh"
echo "Started calculating heterozygous indvidual numbers..."
excavate_output_dir=$output_dir"/excavate/excavate_outputs"
filtered_vcf_dir=$output_dir"/excavate/Guide_filtered_vcfs"
qsub -l mem_free=2G -l h_rt=01:00:00 -sync y "$get_targeted_hets" "$output_dir/excavate/het_individuals" "$param_file" "$genes_w_guides" "$excavate_output_dir" "$filtered_vcf_dir"
echo "Finished calculating heterozygous individual numbers."

# script to calculate number of guides needed to reach heterozygous individuals for each gene
get_guide_info="$project_root/scripts/get_guides/non_excision_guides.sh"
echo "Calculating number of guides to target heterozygotes..."
excavate_output_dir=$output_dir"/excavate/excavate_outputs"
filtered_vcf_dir=$output_dir"/excavate/Guide_filtered_vcfs"
qsub -l mem_free=2G -l h_rt=01:00:00 -sync y "$get_guide_info" "$output_dir/excavate/guide_numbers" "$param_file" "$genes_w_guides" "$excavate_output_dir" "$filtered_vcf_dir"
echo "Finished calculating number of guides."
