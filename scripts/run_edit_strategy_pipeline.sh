#!/bin/bash
#$ -N run_edit_strategy_pipeline
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o logs/out/run_edit_strategy_pipeline.out
#$ -e logs/err/run_edit_strategy_pipeline.err

set -euo pipefail

# Resolve paths relative to this script instead of using hard-coded absolute paths.
project_root="/wynton/home/capra/gramey02/dnd_project"
script_dir="$project_root/scripts"

# source params for the run
param_file="$project_root/data/params/params.txt"
source "$param_file"
echo "$RUN_NAME"

# ---- helper function: create a directory if missing, print a friendly message
mkdir_if_missing() {
  mkdir -p "$1"
}

# # make output directory for the run
output_dir="${OUTPUT_DIR}${RUN_NAME}"
mkdir_if_missing "$output_dir"
param_dir=$output_dir"/PARAMS"
mkdir_if_missing "$param_dir"
cp "$param_file" "$param_dir" # save param file to the output directory for the run (so we can check later if needed)
new_param_file=$output_dir"/PARAMS/params.txt"

# run script to create the necessary directories for the results (right now, set up to have a directory per editing strategy)
make_dirs_script="$script_dir/structure_result_dir/make_dirs.sh"
bash "$make_dirs_script" "$new_param_file"

# run transcript filtering if param is set to true
if (( $FILTER_TRANSCRIPTS == 1 )); then
    mkdir -p "$output_dir/filtered_transcripts" # mkdir -p "$output_dir/filtered_transcripts/expression_filtered"
    # get the number of genes of interest for this run
    unique_genes_file="$output_dir/unique_gene_in_run.txt"
    awk -F',' 'NR > 1 {gsub(/"/, "", $1); gsub(/"/, "", $3); print $1 "\t" $3}' $ORIGINAL_EXON_FILE | sort -u > $unique_genes_file
    # get the number of unique genes in the file
    num_unique_genes=$(wc -l < "$unique_genes_file")
    echo "Number of unique gene/Ensembl ID combinations: $num_unique_genes"
    # run transcript filtering
    if [ "$num_unique_genes" -gt 1000 ]; then
        # since the number of genes is so large, run a parallelized script for transcript filtering
        mkdir -p "$output_dir/filtered_transcripts/per_gene_files"
        filter_script="$script_dir/filter_transcripts/filter_transcripts_all_pc_genes.sh"
        qsub -t 1-$num_unique_genes -l mem_free=20G -l h_rt=03:00:00 -sync y -o "$project_root/logs/out/filter_transcripts.out" -e "$project_root/logs/err/filter_transcripts.err" "$filter_script" "$output_dir" "$new_param_file" "$unique_genes_file"

        # run a short code block to merge all of the information back into one file
        final_filtered_file="$output_dir/filtered_transcripts/filtered_exon_info.csv"
        tmp_body="$output_dir/filtered_transcripts/filtered_exon_info.body.tmp"
        tmp_header="$output_dir/filtered_transcripts/filtered_exon_info.header.tmp"
        sort_tmp="$output_dir/filtered_transcripts/sort_tmp"
        mkdir -p "$sort_tmp"
        # Save the header from the first file
        head -n 1 "$(ls "$output_dir"/filtered_transcripts/per_gene_files/*_filtered_exon_info.csv | head -n 1)" > "$tmp_header"
        # Merge files, skip headers, sort/deduplicate using disk-backed sort
        awk 'FNR==1 {next} {print}' \
            "$output_dir"/filtered_transcripts/per_gene_files/*_filtered_exon_info.csv \
            | sort -u -T "$sort_tmp" -S 4G \
            > "$tmp_body"
        # Reattach header
        cat "$tmp_header" "$tmp_body" > "$final_filtered_file"
        # Clean up
        rm "$tmp_header" "$tmp_body"
        echo "Merged file written to: $final_filtered_file"
        # append to parameter file
        printf 'EXON_FILE_FOR_ANALYSIS="%s"\n' "$final_filtered_file" >> "$new_param_file" # append KEY=VALUE
    else
        filter_script="$script_dir/filter_transcripts/filter_transcripts.sh" # include this in the params file eventually
        qsub -l mem_free=30G -l h_rt=03:00:00 -sync y -o "$project_root/logs/out/filter_transcripts.out" -e "$project_root/logs/err/filter_transcripts.err" "$filter_script" "$output_dir" "$new_param_file"
    fi
else
    # update
    printf 'EXON_FILE_FOR_ANALYSIS="%s"\n' "$ORIGINAL_EXON_FILE" >> "$new_param_file" # append KEY=VALUE
fi

# uncomment below

# now call shell scripts to run sub-pipelines for each of the editing strategies and pass in the set parameters
# parallelized so they can run at the same time

# helper function to parse edit strategies array
contains() {
    local seeking="$1"; shift
    for element in "$@"; do
        [[ "$element" == "$seeking" ]] && return 0
    done
    return 1
}

# indel pipeline
if contains "indels" "${EDIT_STRATS[@]}"; then
    echo "Running indel pipeline..."
    indel_pipeline="$script_dir/edit_strategy_pipelines/indel_pipeline.sh"
    indel_output_dir=$output_dir"/indels"
    qsub -cwd -l mem_free=1G -l h_rt=12:00:00 -o "$project_root/logs/out/indel_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/indel_pipeline_${RUN_NAME}.err" "$indel_pipeline" "$indel_output_dir" "$new_param_file"
    echo "Finished running indel pipeline."
fi

# crisproff pipeline
if contains "CRISPRoff" "${EDIT_STRATS[@]}"; then
    echo "Running crisproff pipeline..."
    crisproff_pipeline="$script_dir/edit_strategy_pipelines/crisproff_pipeline.sh"
    crisproff_output_dir=$output_dir"/CRISPRoff"
    qsub -cwd -l mem_free=1G -l h_rt=12:00:00 -o "$project_root/logs/out/crisproff_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/crisproff_pipeline_${RUN_NAME}.err" "$crisproff_pipeline" "$crisproff_output_dir" "$new_param_file"
    echo "Finished running crisproff pipeline."
fi

# acceptor base edits pipeline
if contains "acceptor_base_edits" "${EDIT_STRATS[@]}"; then
    echo "Running acceptor base edits pipeline..."
    acceptor_pipeline="$script_dir/edit_strategy_pipelines/acceptor_baseEdit_pipeline.sh"
    acceptor_output_dir=$output_dir"/acceptor_base_edits"
    qsub -cwd -l mem_free=1G -l h_rt=16:00:00 -o "$project_root/logs/out/acceptor_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/acceptor_pipeline_${RUN_NAME}.err" "$acceptor_pipeline" "$acceptor_output_dir" "$new_param_file"
    echo "Finished running acceptor pipeline."
fi

# donor base edits pipeline
if contains "donor_base_edits" "${EDIT_STRATS[@]}"; then
    echo "Running donor base edits pipeline..."
    donor_pipeline="$script_dir/edit_strategy_pipelines/donor_baseEdit_pipeline.sh"
    donor_output_dir=$output_dir"/donor_base_edits"
    qsub -cwd -l mem_free=1G -l h_rt=12:00:00 -o "$project_root/logs/out/donor_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/donor_pipeline_${RUN_NAME}.err" "$donor_pipeline" "$donor_output_dir" "$new_param_file"
    echo "Finished running donor pipeline."
fi

# excision pipeline
if contains "excision" "${EDIT_STRATS[@]}"; then
    echo "Running excision pipeline..."
    excision_pipeline="$script_dir/edit_strategy_pipelines/excision_pipeline.sh"
    excision_output_dir=$output_dir"/excision"
    qsub -cwd -l mem_free=1G -l h_rt=16:00:00 -o "$project_root/logs/out/excision_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/excision_pipeline_${RUN_NAME}.err" "$excision_pipeline" "$excision_output_dir" "$new_param_file"
    echo "Finished running excision pipeline."
fi
