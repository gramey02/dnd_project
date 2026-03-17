#!/bin/bash
#$ -N make_dirs
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/protected/home/capra/gramey02/ConklinCollab/scripts/out/make_dirs.out
#$ -e /wynton/protected/home/capra/gramey02/ConklinCollab/scripts/err/make_dirs.err

# script to create directories for the editing pipeline results

# load input argument that contains parameters for the run
param_file=$1
source $param_file
# overall output directory for the run
output_dir=$OUTPUT_DIR$RUN_NAME

# ---- helper function: create a directory if missing, print a friendly message
mkdir_if_missing() {
  local d="$1"
  if [ -d "$d" ]; then
    echo "Directory '$d' already exists."
  else
    mkdir -p "$d"       # -p creates parents as needed; safe to call repeatedly
    echo "Directory '$d' created."
  fi
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















# --------------------------------Old code, made more efficient above--------------------

# --------------------------------------------------------------indels-------------------------------------------
# if [[ " ${EDIT_STRATS[@]} " =~ " indels " ]]; then
#     # overall directory for the specific editing strategy
#     DIRECTORY=$output_dir"/indels"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for ubiquitous regions
#     DIRECTORY=$output_dir"/indels/ubiq_regions"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for common vars in those regions
#     DIRECTORY=$output_dir"/indels/ubiq_region_CommonVars"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for NMD predictions within exons (for this specific editing strategy)
#     DIRECTORY=$output_dir"/indels/NMD"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for guide-filtered snp positions
#     DIRECTORY=$output_dir"/indels/excavate/Guide_locs"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for guide-filtered snp positions
#     DIRECTORY=$output_dir"/indels/excavate/Guide_filtered_vcfs"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for heterozygous individual calculations
#     DIRECTORY=$output_dir"/indels/excavate/het_individuals"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for heterozygous individual calculations
#     DIRECTORY=$output_dir"/indels/excavate/het_individuals/metadata"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # overall directory for excavate readouts
#     DIRECTORY=$output_dir"/indels/excavate"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for specific excavate readouts/inputs
#     DIRECTORY=$output_dir"/indels/excavate/CommonVar_locs"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi
#     DIRECTORY=$output_dir"/indels/excavate/input_metadata"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi
#     DIRECTORY=$output_dir"/indels/excavate/input_vcfs"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi
#     DIRECTORY=$output_dir"/indels/excavate/excavate_outputs"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi
# fi

# # --------------------------------------------------------base editing--------------------------------------------------
# if [[ " ${EDIT_STRATS[@]} " =~ " base_edits " ]]; then
#     # overall directory for the specific editing strategy
#     DIRECTORY=$output_dir"/base_edits"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for ubiquitous regions
#     DIRECTORY=$output_dir"/base_edits/donor_sites/ubiq_regions"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for common vars in those regions
#     DIRECTORY=$output_dir"/base_edits/donor_sites/ubiq_region_CommonVars"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for excavate readouts
#     DIRECTORY=$output_dir"/base_edits/donor_sites/excavate"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # do the same below for acceptor regions
#     DIRECTORY=$output_dir"/base_edits/acceptor_sites/ubiq_regions"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi
#     DIRECTORY=$output_dir"/base_edits/acceptor_sites/ubiq_region_CommonVars"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

# fi

# # ------------------------------------CRISPRoff--------------------------------------
# if [[ " ${EDIT_STRATS[@]} " =~ " CRISPRoff " ]]; then
#     # overall directory for the specific editing strategy
#     DIRECTORY=$output_dir"/CRISPRoff"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for ubiquitous regions
#     DIRECTORY=$output_dir"/CRISPRoff/ubiq_regions"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for common vars in those regions
#     DIRECTORY=$output_dir"/CRISPRoff/ubiq_region_CommonVars"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

#     # directory for excavate readouts
#     DIRECTORY=$output_dir"/CRISPRoff/excavate"
#     # Check if the directory exists using the -d test operator; if not, create it
#     if [ ! -d "$DIRECTORY" ]; then
#         mkdir -p "$DIRECTORY"
#         echo "Directory '$DIRECTORY' created."
#     else
#         echo "Directory '$DIRECTORY' already exists."
#     fi

# fi

# # ----------------------------------------excision-----------------------------------------
# # -----------------------------------Chromosome Set Directory------------------------------
# DIRECTORY=$output_dir"/chromosomes"
# # Check if the directory exists using the -d test operator; if not, create it
# if [ ! -d "$DIRECTORY" ]; then
#     mkdir -p "$DIRECTORY"
#     echo "Directory '$DIRECTORY' created."
# else
#     echo "Directory '$DIRECTORY' already exists."
# fi
# # -----------------------------------Summary Files Directory-------------------------------
# DIRECTORY=$output_dir"/summary_files"
# # Check if the directory exists using the -d test operator; if not, create it
# if [ ! -d "$DIRECTORY" ]; then
#     mkdir -p "$DIRECTORY"
#     echo "Directory '$DIRECTORY' created."
# else
#     echo "Directory '$DIRECTORY' already exists."
# fi