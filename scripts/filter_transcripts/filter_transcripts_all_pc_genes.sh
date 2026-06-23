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
unique_genes="$3"

# run script below

# isolate each gene
ensembl_id=$(awk -v row="$SGE_TASK_ID" 'NR == row {id=$1; gsub(/"/,"",id); sub(/\.[0-9]+$/,"",id); print id}' "$unique_genes")
gene=$(awk -v row="$SGE_TASK_ID" 'NR == row {g=$2; gsub(/"/,"",g); print g}' "$unique_genes")

# isolate the transcript expression for this gene only
mkdir -p "$project_root/data/dnd_ensembl/per_gene_files"
cur_transcript_tpm_file="$project_root/data/dnd_ensembl/per_gene_files/${gene}_${ensembl_id}_expression.txt"
awk -v target_id="$ensembl_id" '{cur_id=$2; gsub(/"/,"",cur_id); sub(/\.[0-9]+$/,"",cur_id); if (cur_id == target_id) print}' "$project_root/$TRANSCRIPT_TPM_FILE" > "$cur_transcript_tpm_file"

# pass the filtered expression data into the python script
script="$script_dir/filter_transcripts/filter_transcripts_expression_all_pc_genes.py"
python3 $script --transcript_tpm_file $cur_transcript_tpm_file \
    --sample_attributes_file "$project_root/$SAMPLE_ATTRIBUTES_FILE" \
    --gene_median_tpms_file "$project_root/$GENE_MEDIAN_TPMS_FILE" \
    --exon_file "$project_root/$cur_exon_file" \
    --tpm_thresh "$GENE_EXPRESSION_THRESH" \
    --prop_thresh "$EXPRESSION_PROP" \
    --keep_all_transcripts "$KEEP_ALL_TRANSCRIPTS" \
    --output_file "$output_dir/filtered_transcripts/per_gene_files/${gene}_${ensembl_id}_filtered_exon_info.csv" \
    --colnames_file "$project_root/$COLNAMES_FILE" \
    --tissue_map_file "$project_root/$TISSUE_MAP_FILE"