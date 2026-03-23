#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

# parse input arguments
output_dir="$1"
param_file="$2"
source "$param_file"
cur_exon_file="$ORIGINAL_EXON_FILE"

# Resolve repo-relative paths from params.txt against the project root.
if [[ "$cur_exon_file" != /* ]]; then
    cur_exon_file="$(cd "$project_root" && realpath "$cur_exon_file")"
fi

# run script below
echo "Filtering transcripts by expression proportion..."
# run the expression-based transcript filtering script
expr_filt_script="$script_dir/filter_transcripts_expression.py"
python3 "$expr_filt_script" --transcript_tpm_file "$TRANSCRIPT_TPM_FILE" \
    --sample_attributes_file "$SAMPLE_ATTRIBUTES_FILE" \
    --gene_median_tpms_file "$GENE_MEDIAN_TPMS_FILE" \
    --exon_file "$cur_exon_file" \
    --tpm_thresh "$GENE_EXPRESSION_THRESH" \
    --prop_thresh "$EXPRESSION_PROP" \
    --keep_all_transcripts "$KEEP_ALL_TRANSCRIPTS" \
    --output_file "$output_dir/filtered_transcripts/filtered_exon_info.csv"
echo "Finished filtering transcripts by expression proportion."
cur_exon_file="$output_dir/filtered_transcripts/filtered_exon_info.csv"

# # below, filter on tsl or canonical flags if true
# if [[ $CANONICAL_TSL_FILTER == "True" ]]; then
#     echo "Filtering transcripts by tsl and canonical flags..."
#     canonical_tsl_filt_script="/wynton/protected/home/capra/gramey02/ConklinCollab/scripts/DN_CommonVars/filter_transcripts/filter_transcripts_tsl_canonical.py"
#     python3 $canonical_tsl_filt_script --use_only_canonical $USE_ONLY_CANONICAL \
#         --use_tsl_flag $USE_TSL_FLAG \
#         --tsl_num $TSL_NUM \
#         --exon_file $cur_exon_file \
#         --output_file "$ouptut_dir/filtered_transcripts/tsl_canonical_filtered/exon_info.csv"
#     cur_exon_file="$ouptut_dir/filtered_transcripts/tsl_canonical_filtered/exon_info.csv"
# fi

# add final filtered exon file to params file for future use
final_exon_file="$cur_exon_file"
printf 'EXON_FILE_FOR_ANALYSIS="%s"\n' "$final_exon_file" >> "$param_file" # append KEY=VALUE
