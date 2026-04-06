#!/bin/bash
#$ -N filter_transcripts
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input arguments
output_dir="$1"
param_file="$2"
source "$param_file"
cur_exon_file="$ORIGINAL_EXON_FILE"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# run script below
echo "Filtering transcripts by expression proportion..."
# run the expression-based transcript filtering script
expr_filt_script="$script_dir/filter_transcripts/filter_transcripts_expression.py"
python3 "$expr_filt_script" --transcript_tpm_file "$project_root/$TRANSCRIPT_TPM_FILE" \
    --sample_attributes_file "$project_root/$SAMPLE_ATTRIBUTES_FILE" \
    --gene_median_tpms_file "$project_root/$GENE_MEDIAN_TPMS_FILE" \
    --exon_file "$project_root/$cur_exon_file" \
    --tpm_thresh "$GENE_EXPRESSION_THRESH" \
    --prop_thresh "$EXPRESSION_PROP" \
    --keep_all_transcripts "$KEEP_ALL_TRANSCRIPTS" \
    --output_file "$output_dir/filtered_transcripts/filtered_exon_info.csv" \
    --colnames_file "$project_root/$COLNAMES_FILE" \
    --tissue_map_file "$project_root/$TISSUE_MAP_FILE"
echo "Finished filtering transcripts by expression proportion."
cur_exon_file="$output_dir/filtered_transcripts/filtered_exon_info.csv"

# add final filtered exon file to params file for future use
final_exon_file="$cur_exon_file"
printf 'EXON_FILE_FOR_ANALYSIS="%s"\n' "$final_exon_file" >> "$param_file" # append KEY=VALUE
