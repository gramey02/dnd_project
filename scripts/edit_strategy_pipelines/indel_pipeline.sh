#!/bin/bash
#$ -N indel_pipeline
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/indel_pipeline.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/indel_pipeline.err

# parse input arguments
param_file=$2
source $param_file
output_dir=$1

# find exon_file
#exon_file=$EXON_FILE_FOR_ANALYSIS
#exon_file="/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/filtered_exon_info/test_exon.csv"
#exon_file="/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/filtered_exon_info/dHS_exonInfo_GTExCurated_transcripts_TPM1.0_ExprProp0.01_2025_04_22.csv"
exon_file=$EXON_FILE_FOR_ANALYSIS

# get num chromosomes
col=$(head -1 $exon_file | tr ',' '\n' | grep -nx "chromosome_name" | cut -d: -f1)
cut -d',' -f"$col" $exon_file | tail -n +2 | sort -u > $OUTPUT_DIR"/"$RUN_NAME"/chromosomes/chrom_set.txt"
num_chroms=$(wc -l < $OUTPUT_DIR"/"$RUN_NAME"/chromosomes/chrom_set.txt")

# script to get ubiquitous exonic regions & identify common vars in them
ubiq_regions_CommonVar_script="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_CommonVars/get_common_vars_UER.sh"
echo "Started identifying common coding vars..."
qsub -l mem_free=1G -l h_rt=00:20:00 -t 1-"$num_chroms" -sync y $ubiq_regions_CommonVar_script "$output_dir" "$param_file" "$num_chroms" "$exon_file"
echo "Finished identifiying common coding vars."

##### NMD script goes here
targetable_common_var_dict="$output_dir/ubiq_region_CommonVars/CommonVars_ALL_dict.pkl"
NMD_script="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_NMD/Annot_NMD_escape.sh"
echo "Started running NMD analysis on guide-filtered indel vars..."
qsub -l mem_free=1G -l h_rt=01:00:00 -sync y $NMD_script "$output_dir" "$param_file" "$targetable_common_var_dict"
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
excavate_vcf_creation="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_PAMsites/generate_filtered_vcfs.sh" # first we need to generate vcf.gz files for the genes that have variants in their ubiquitous exonic regions
input_metadata=$output_dir"/excavate/input_metadata/excavate_run_metadata.txt"
echo "Started creating vcf files for excavate input..."
qsub -t 1-"$num_common_var_genes" -l mem_free=2G -l h_rt=01:00:00 -sync y $excavate_vcf_creation "$output_dir" "$param_file" "$input_metadata"
echo "Finished creating excavate inputs."

# script to run excavate
run_excavate_script="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_PAMsites/excavate/run_excavate.sh"
echo "Started running EXCAVATE..."
qsub -t 1-"$num_common_var_genes" -l mem_free=1G -l h_rt=00:45:00 -sync y $run_excavate_script "$output_dir" "$param_file" "$input_metadata"
echo "Finished running EXCAVATE."

# generate text files for the valid guides so you can filter the vcfs accordingly
generate_guide_textFiles="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_PAMsites/generate_guide_textFiles.py"
guides_filepath=$output_dir"/excavate/excavate_outputs" # switch to not include indels later
echo "Started generating position files for valid guides..."
python3 $generate_guide_textFiles --guides_filepath $guides_filepath --exon_file $exon_file --output_dir $output_dir
echo "Finished generating position files for valid guides."

# filter vcfs based on these viable EXCAVATE guides to see how many heterozygous individuals they will capture
guide_based_filtering="/wynton/home/capra/gramey02/ConklinCollab/scripts/bcftools_examples/position_filtering.sh"
genes_w_guides=$output_dir"/excavate/het_individuals/metadata/genes_w_valid_guides.txt"
num_genes_w_guides=$(awk -F'\t' '$1 != "" {n++} END{print n}' $genes_w_guides)
echo "Started filtering vcfs based on valid guides..."
qsub -t 1-"$num_genes_w_guides" -l mem_free=2G -l h_rt=01:00:00 -sync y $guide_based_filtering "$output_dir" "$param_file" "$genes_w_guides"
echo "Finished filtering vcfs based on valid guides."

# script to calculate number of heterozygous individuals for each gene
get_targeted_hets="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_targeted_hets/get_targeted_hets.sh"
echo "Started calculating heterozygous indvidual numbers..."
excavate_output_dir=$output_dir"/excavate/excavate_outputs"
filtered_vcf_dir=$output_dir"/excavate/Guide_filtered_vcfs"
qsub -l mem_free=2G -l h_rt=01:00:00 -sync y $get_targeted_hets "$output_dir/excavate/het_individuals" "$param_file" "$genes_w_guides" "$excavate_output_dir" "$filtered_vcf_dir"
echo "Finished calculating heterozygous individual numbers."

# script to calculate number of guides needed to reach heterozygous individuals for each gene
get_guide_info="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/get_guides.sh"
echo "Calculating number of guides to target heterozygotes..."
excavate_output_dir=$output_dir"/excavate/excavate_outputs"
filtered_vcf_dir=$output_dir"/excavate/Guide_filtered_vcfs"
qsub -l mem_free=2G -l h_rt=01:00:00 -sync y $get_guide_info "$output_dir/excavate/guide_numbers" "$param_file" "$genes_w_guides" "$excavate_output_dir" "$filtered_vcf_dir"
echo "Finished calculating number of guides."






#-----------------code to generate multiple sets of variants----------------
# # find common vars slightly differently if we're using shared regions among transcripts or if we want to consider all common vars. Store each created set in an array
# variant_sets=()
# summary_file_sets=()
# branch_names=()
# sub_branch_names=()

# if (( $USE_ALL_REGIONS == 1 )); then
#     # script to get common vars across transcripts, no matter the ubiquitous regions
#     CommonVar_script="/wynton/protected/home/capra/gramey02/ConklinCollab/scripts/DN_CommonVars/get_common_vars_tp.sh"
#     echo "Started identifying common coding vars across all transcript exons..."
#     qsub -l mem_free=1G -l h_rt=00:30:00 -t 1-"$num_chroms" -sync y $CommonVar_script "$output_dir/all_regions" "$param_file" "$num_chroms" "$exon_file"
#     echo "Finished identifying common coding vars across transcripts."

#     # add to variant set array
#     # need to change these to have top level file name in there
#     variant_sets+=( "$output_dir/all_regions/all_region_CommonVars/CommonVars_ALL_dict.pkl" )
#     summary_file_sets+=( "$output_dir/all_regions/all_region_CommonVars/CommonVars_ALL_summary_noIDX.txt" )
#     branch_names+=( "$output_dir/all_regions" )
#     sub_branch_names+=( "$output_dir/all_regions/all_region_CommonVars")
# fi

# if (( $USE_UBIQUITOUS_REGIONS == 1 )); then
#     # script to get ubiquitous exonic regions & identify common vars in them
#     ubiq_regions_CommonVar_script="/wynton/protected/home/capra/gramey02/ConklinCollab/scripts/DN_CommonVars/get_common_vars_UER.sh" # UER = ubiquitous exonic regions
#     echo "Started identifying common coding vars in ubiquitous exonic regions..."
#     qsub -l mem_free=1G -l h_rt=00:20:00 -t 1-"$num_chroms" -sync y $ubiq_regions_CommonVar_script "$output_dir/ubiq_regions" "$param_file" "$num_chroms" "$exon_file"
#     echo "Finished identifying common coding vars in ubiquitous regions."

#     # add to variant set array
#     # need to change these to have top level file name in thre (ubiq_regions/ubiq_region_CommonVars)
#     variant_sets+=( "$output_dir/ubiq_regions/ubiq_region_CommonVars/CommonVars_ALL_dict.pkl" )
#     summary_file_sets+=( "$output_dir/ubiq_regions/ubiq_region_CommonVars/CommonVars_ALL_summary_noIDX.txt" )
#     branch_names+=( "$output_dir/ubiq_regions" )
#     sub_branch_names+=( "$output_dir/ubiq_regions/ubiq_region_CommonVars")

# fi


# function runs downstream analyses

# for i in "${!variant_sets[@]}"; do
#     variant_set="${variant_sets[$i]}"
#     summary_file="${summary_file_sets[$i]}"
#     branch_name="${branch_names[$i]}"
#     sub_branch_name="${sub_branch_names[$i]}"

#     # script to run downstream analysis
#     downstream_analysis_script="/wynton/protected/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/run_downstream_analysis_indels.sh"
    
#     # input and output files
#     var_set_identifier=$(basename "$branch_name")
#     job_output="/wynton/protected/home/capra/gramey02/ConklinCollab/scripts/out/run_downstream_analysis_indels_${var_set_identifier}.out"
#     job_err="/wynton/protected/home/capra/gramey02/ConklinCollab/scripts/out/run_downstream_analysis_indels_${var_set_identifier}.err"

#     qsub -l mem_free=1G -l h_rt=03:00:00 \
#         -o $job_output \
#         -e $job_err \
#         $downstream_analysis_script \
#         "$variant_set" "$summary_file" "$branch_name" "$sub_branch_name" "$param_file"
# done 

## Old NMD code
# # script to run NMD pipeline on these guide-filtered common vars
# NMD_script="/wynton/protected/home/capra/gramey02/ConklinCollab/scripts/DN_NMD/Annot_NMD_escape.sh"
# targetable_common_var_dict="$output_dir/excavate/het_individuals/metadata/targetable_SNP_locs_by_gene.pkl"
# echo "Started running NMD analysis on guide-filtered indel vars..."
# qsub -l mem_free=1G -l h_rt=01:00:00 -sync y $NMD_script "$output_dir/NMD" "$param_file" "$targetable_common_var_dict"
# echo "Finished running NMD analysis."