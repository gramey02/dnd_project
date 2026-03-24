#!/bin/bash
#$ -N make_dirs
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ../../logs/out/make_dirs.out
#$ -e ../../logs/err/make_dirs.err

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

# script to create directories for the editing pipeline results
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"
execution_utils="$project_root/scripts/utils/execution_mode.sh"

# load input argument that contains parameters for the run
param_file="$1"
source "$param_file"
source "$execution_utils"

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

# Create HPC log directories at the run level when requested.
run_mode="$(get_run_mode)"
if [[ "$run_mode" == "hpc" ]]; then
  mkdir_if_missing "$output_dir/logs"
  mkdir_if_missing "$output_dir/logs/out"
  mkdir_if_missing "$output_dir/logs/err"
fi

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
