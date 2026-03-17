#!/bin/bash
#$ -N get_targeted_hets_prePAM
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/get_targeted_hets_prePAM.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/get_targeted_hets_prePAM.err

# parse input arguments
param_file="/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/params.txt"
source $param_file

mkdir_if_missing() {
  local d="$1"
  if [ -d "$d" ]; then
    echo "Directory '$d' already exists."
  else
    mkdir -p "$d"       # -p creates parents as needed; safe to call repeatedly
    echo "Directory '$d' created."
  fi
}

# # indel inputs
# cur_strat="indels"
# # # crisproff inputs
# # cur_strat="CRISPRoff"
# # # donor_baseEdit inputs
# # cur_strat="donor_base_edits"
# # # acceptor_baseEdit inputs
# # cur_strat="acceptor_base_edits"

# # set up directories and variables
# make_het_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/prePAM_hets"
# mkdir_if_missing $make_het_dir
# output_dir=$make_het_dir
# common_var_vcf_file_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_vcfs"
# gene_info="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_metadata/excavate_run_metadata.txt"

# echo $make_het_dir

# # unzip files for input
# for f in "$common_var_vcf_file_dir"/*.vcf.gz; do
#     out="${f%.gz}"

#     if [[ ! -e "$out" ]]; then
#         gzip -dk "$f"
#     fi
# done

# # # run het individual calculation script
# # # (script will save both the number of hets and the identifiers of the hets
# # # to be compared to other editing methods)
# script="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_targeted_hets/get_targeted_hets_prePAM.py"

# python3 $script --output_dir $output_dir \
#     --gene_info $gene_info \
#     --filtered_vcf_dir $common_var_vcf_file_dir \
#     --num_samples $NUM_SAMPLES


# # ------------------ crisproff ------------------- #
# # crisproff inputs
# cur_strat="CRISPRoff"
# # # donor_baseEdit inputs
# # cur_strat="donor_base_edits"
# # # acceptor_baseEdit inputs
# # cur_strat="acceptor_base_edits"

# # set up directories and variables
# make_het_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/prePAM_hets"
# mkdir_if_missing $make_het_dir
# output_dir=$make_het_dir
# common_var_vcf_file_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_vcfs"
# gene_info="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_metadata/excavate_run_metadata.txt"

# echo $make_het_dir

# # unzip files for input
# for f in "$common_var_vcf_file_dir"/*.vcf.gz; do
#     out="${f%.gz}"

#     if [[ ! -e "$out" ]]; then
#         gzip -dk "$f"
#     fi
# done

# # # run het individual calculation script
# # # (script will save both the number of hets and the identifiers of the hets
# # # to be compared to other editing methods)
# script="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_targeted_hets/get_targeted_hets_prePAM.py"

# python3 $script --output_dir $output_dir \
#     --gene_info $gene_info \
#     --filtered_vcf_dir $common_var_vcf_file_dir \
#     --num_samples $NUM_SAMPLES

# # ----------- donor base edits ------------------------
# # donor_baseEdit inputs
# cur_strat="donor_base_edits"
# # # acceptor_baseEdit inputs
# # cur_strat="acceptor_base_edits"

# # set up directories and variables
# make_het_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/prePAM_hets"
# mkdir_if_missing $make_het_dir
# output_dir=$make_het_dir
# common_var_vcf_file_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_vcfs"
# gene_info="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_metadata/excavate_run_metadata.txt"

# echo $make_het_dir

# # unzip files for input
# for f in "$common_var_vcf_file_dir"/*.vcf.gz; do
#     out="${f%.gz}"

#     if [[ ! -e "$out" ]]; then
#         gzip -dk "$f"
#     fi
# done

# # # run het individual calculation script
# # # (script will save both the number of hets and the identifiers of the hets
# # # to be compared to other editing methods)
# script="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_targeted_hets/get_targeted_hets_prePAM.py"

# python3 $script --output_dir $output_dir \
#     --gene_info $gene_info \
#     --filtered_vcf_dir $common_var_vcf_file_dir \
#     --num_samples $NUM_SAMPLES


# # ---------------------- acceptor base edits ---------------
# # acceptor_baseEdit inputs
# cur_strat="acceptor_base_edits"

# # set up directories and variables
# make_het_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/prePAM_hets"
# mkdir_if_missing $make_het_dir
# output_dir=$make_het_dir
# common_var_vcf_file_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_vcfs"
# gene_info="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_metadata/excavate_run_metadata.txt"

# echo $make_het_dir

# # unzip files for input
# for f in "$common_var_vcf_file_dir"/*.vcf.gz; do
#     out="${f%.gz}"

#     if [[ ! -e "$out" ]]; then
#         gzip -dk "$f"
#     fi
# done

# # # run het individual calculation script
# # # (script will save both the number of hets and the identifiers of the hets
# # # to be compared to other editing methods)
# script="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_targeted_hets/get_targeted_hets_prePAM.py"

# python3 $script --output_dir $output_dir \
#     --gene_info $gene_info \
#     --filtered_vcf_dir $common_var_vcf_file_dir \
#     --num_samples $NUM_SAMPLES


strategies=("indels" "CRISPRoff" "donor_base_edits" "acceptor_base_edits")

for cur_strat in "${strategies[@]}"; do
    make_het_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/prePAM_hets"
    mkdir_if_missing "$make_het_dir"
    output_dir=$make_het_dir
    common_var_vcf_file_dir="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_vcfs"
    gene_info="$OUTPUT_DIR$RUN_NAME/$cur_strat/excavate/input_metadata/excavate_run_metadata.txt"

    echo "Processing $cur_strat..."
    echo "Output dir: $make_het_dir"

    # unzip VCFs if needed
    for f in "$common_var_vcf_file_dir"/*.vcf.gz; do
        out="${f%.gz}"
        [[ -e "$out" ]] && continue
        gzip -dk "$f" &>/dev/null
    done

    python3 /wynton/home/capra/gramey02/ConklinCollab/scripts/DN_targeted_hets/get_targeted_hets_prePAM.py \
        --output_dir "$output_dir" \
        --gene_info "$gene_info" \
        --filtered_vcf_dir "$common_var_vcf_file_dir" \
        --num_samples "$NUM_SAMPLES"
done