#!/bin/bash
#$ -N get_targeted_hets_prePAM
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input arguments
output_dir="$1"
param_file="$2"
gene_info="$3"
filtered_vcf_dir="$4"
source "$param_file"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# ensure the output dir exists
mkdir -p $output_dir

# run script
script="$script_dir/get_hets/get_targeted_hets_prePAM.py"
python3 "$script" --output_dir "$output_dir" --gene_info "$gene_info" --filtered_vcf_dir "$filtered_vcf_dir" --num_samples "$NUM_SAMPLES"











# strategies=("indels" "CRISPRoff" "donor_base_edits" "acceptor_base_edits")

# for cur_strat in "${strategies[@]}"; do
#     make_het_dir="$resolved_output_base$RUN_NAME/$cur_strat/prePAM_hets"
#     mkdir_if_missing "$make_het_dir"
#     output_dir="$make_het_dir"
#     common_var_vcf_file_dir="$resolved_output_base$RUN_NAME/$cur_strat/excavate/input_vcfs"
#     gene_info="$resolved_output_base$RUN_NAME/$cur_strat/excavate/input_metadata/excavate_run_metadata.txt"

#     echo "Processing $cur_strat..."
#     echo "Output dir: $make_het_dir"

#     # unzip VCFs if needed
#     for f in "$common_var_vcf_file_dir"/*.vcf.gz; do
#         out="${f%.gz}"
#         [[ -e "$out" ]] && continue
#         gzip -dk "$f" &>/dev/null
#     done

#     python3 "$script_dir/get_targeted_hets_prePAM.py" \
#         --output_dir "$output_dir" \
#         --gene_info "$gene_info" \
#         --filtered_vcf_dir "$common_var_vcf_file_dir" \
#         --num_samples "$NUM_SAMPLES"
