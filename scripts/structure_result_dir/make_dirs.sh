#!/bin/bash
#$ -N make_dirs
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ../../logs/out/make_dirs.out
#$ -e ../../logs/err/make_dirs.err

# load input argument that contains parameters for the run
param_file="$1"
source "$param_file"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

# overall output directory for the run
output_dir="${resolved_output_base}${RUN_NAME}"

# ---- helper function: create a directory if missing, print a friendly message
mkdir_if_missing() {
  local d="$1"
  mkdir -p "$d"
}

# ---- 1) Always create these at the RUN LEVEL (outside editing strategies)
ALWAYS_TOPLEVEL=( "summary_files" "chromosomes" )

# Make sure the overall run dir itself exists
mkdir_if_missing "$output_dir"

# Create each top-level "always" folder
for d in "${ALWAYS_TOPLEVEL[@]}"; do
  mkdir_if_missing "$output_dir/$d"
done

# ---- 2) Always make these INSIDE EVERY STRATEGY FOLDER
# These are relative to: $output_dir/<strategy>/
# Keeping "" means: also create the strategy root itself.
ALWAYS_PER_STRAT_DEFAULT=(
  ""                             # creates $output_dir/<strategy> (root)
  "ubiq_regions"                 # creates $output_dir/<strategy>/ubiq_regions
  "ubiq_region_CommonVars"       # creates $output_dir/<strategy>/ubiq_regions
  "excavate"                     # etc.
  "excavate/Guide_locs"
  "excavate/Guide_filtered_vcfs"
  "excavate/het_individuals"
  "excavate/het_individuals/metadata"
  "excavate/CommonVar_locs"
  "excavate/input_metadata"
  "excavate/input_vcfs"
  "excavate/excavate_outputs"
  "excavate/guide_numbers"
)

# since excision won't require all of the same directories, we separate its subdirectories out
ALWAYS_PER_STRAT_EXCISION=(
  ""
  "CommonVars"
  "CommonVars/refined_common_vars"
  "CommonVars/valid_snp_pairs"
  "excavate"
  "excavate/Guide_locs"
  "excavate/Guide_filtered_vcfs"
  "excavate/het_individuals"
  "excavate/het_individuals/metadata"
  "excavate/CommonVar_locs"
  "excavate/input_metadata"
  "excavate/input_vcfs"
  "excavate/excavate_outputs"
  "excavate/guide_numbers"
  "final_sgRNA_snps"
  "het_individuals"
)

# ---- 3) Strategy-specific add-ons (only for that strategy)
# All paths here are relative to $output_dir/<strategy>/
declare -a DIRS_indels=(
  "NMD"
)

declare -a DIRS_donor_base_edits=(
    # nothing for now, add too later
)

declare -a DIRS_acceptor_base_edits=(
  # nothing for now, add too later
)

declare -a DIRS_CRISPRoff=(
  "CpG_islands"
  "GC_content"
)

declare -a DIRS_excision=(
  # nothing for now, add to later
)

# ---- 4) For each requested strategy, compose:
#     dirs_to_create = ALWAYS_PER_STRAT + that strategy’s extras
for strat in "${EDIT_STRATS[@]}"; do
  base="$output_dir/$strat"

  # Pick the right "always" set
  case "$strat" in
    excision) always=("${ALWAYS_PER_STRAT_EXCISION[@]}") ;;
    *)         always=("${ALWAYS_PER_STRAT_DEFAULT[@]}") ;;
  esac

  # Create the chosen “always” dirs
  for suffix in "${always[@]}"; do
    mkdir_if_missing "$base${suffix:+/$suffix}"
  done

  # Strategy-specific extras
  case "$strat" in
    indels)                 dirs=("${DIRS_indels[@]}") ;;
    acceptor_base_edits)    dirs=("${DIRS_acceptor_base_edits[@]}") ;;
    donor_base_edits)       dirs=("${DIRS_donor_base_edits[@]}") ;;
    excision)               dirs=("${DIRS_excision[@]}") ;;
    CRISPRoff)              dirs=("${DIRS_CRISPRoff[@]}") ;;
    *)                      dirs=() ;;
  esac
  for suffix in "${dirs[@]}"; do
    mkdir_if_missing "$base/${suffix}"
  done
done
