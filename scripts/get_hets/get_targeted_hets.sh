#!/bin/bash
#$ -N get_targeted_hets
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/get_targeted_hets.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/get_targeted_hets.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

# parse input arguments
output_dir="$1"
param_file="$2"
source "$param_file"
gene_info="$3"
excavate_output_dir="$4"
filtered_vcf_dir="$5"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

# calculate the number of editing strategies you want to combine results for this analysis
strats_to_remove=( "excision" )
declare -A RM=()
for s in "${strats_to_remove[@]}"; do
  RM["${s,,}"]=1
done
filtered=() # Filter EDIT_STRATS -> keep only those NOT in RM (case-insensitive)
for s in "${EDIT_STRATS[@]}"; do
  [[ -n ${RM["${s,,}"]+x} ]] || filtered+=("$s")
done
num_strats=${#filtered[@]} # Count

# run het individual calculation script
# (script will save both the number of hets and the identifiers of the hets
# to be compared to other editing methods)
script="$script_dir/get_targeted_hets.py"
python3 "$script" --output_dir "$output_dir" \
    --gene_info "$gene_info" \
    --excavate_output_dir "$excavate_output_dir" \
    --filtered_vcf_dir "$filtered_vcf_dir" \
    --num_strats "$num_strats" \
    --strats "${filtered[@]}" \
    --run_dir "$resolved_output_base$RUN_NAME"


