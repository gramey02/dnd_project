#!/bin/bash
#$ -N run_edit_strategy_pipeline
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o logs/out/run_edit_strategy_pipeline.out
#$ -e logs/err/run_edit_strategy_pipeline.err

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
    # run transcript filtering
    filter_script="$script_dir/filter_transcripts/filter_transcripts.sh" # include this in the params file eventually
    qsub -l mem_free=1G -l h_rt=01:00:00 -sync y -o "$project_root/logs/out/filter_transcripts.out" -e "$project_root/logs/err/filter_transcripts.err" "$filter_script" "$output_dir" "$new_param_file"
else
    # update
    printf 'EXON_FILE_FOR_ANALYSIS="%s"\n' "$ORIGINAL_EXON_FILE" >> "$new_param_file" # append KEY=VALUE
fi

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
    qsub -cwd -l mem_free=1G -l h_rt=04:00:00 -o "$project_root/logs/out/indel_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/indel_pipeline_${RUN_NAME}.err" "$indel_pipeline" "$indel_output_dir" "$new_param_file"
    echo "Finished running indel pipeline."
fi

# crisproff pipeline
if contains "CRISPRoff" "${EDIT_STRATS[@]}"; then
    echo "Running crisproff pipeline..."
    crisproff_pipeline="$script_dir/edit_strategy_pipelines/crisproff_pipeline.sh"
    crisproff_output_dir=$output_dir"/CRISPRoff"
    qsub -cwd -l mem_free=1G -l h_rt=06:00:00 -o "$project_root/logs/out/crisproff_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/crisproff_pipeline_${RUN_NAME}.err" "$crisproff_pipeline" "$crisproff_output_dir" "$new_param_file"
    echo "Finished running crisproff pipeline."
fi

# acceptor base edits pipeline
if contains "acceptor_base_edits" "${EDIT_STRATS[@]}"; then
    echo "Running acceptor base edits pipeline..."
    acceptor_pipeline="$script_dir/edit_strategy_pipelines/acceptor_baseEdit_pipeline.sh"
    acceptor_output_dir=$output_dir"/acceptor_base_edits"
    qsub -cwd -l mem_free=1G -l h_rt=06:00:00 -o "$project_root/logs/out/acceptor_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/acceptor_pipeline_${RUN_NAME}.err" "$acceptor_pipeline" "$acceptor_output_dir" "$new_param_file"
    echo "Finished running acceptor pipeline."
fi

# donor base edits pipeline
if contains "donor_base_edits" "${EDIT_STRATS[@]}"; then
    echo "Running donor base edits pipeline..."
    donor_pipeline="$script_dir/edit_strategy_pipelines/donor_baseEdit_pipeline.sh"
    donor_output_dir=$output_dir"/donor_base_edits"
    qsub -cwd -l mem_free=1G -l h_rt=06:00:00 -o "$project_root/logs/out/donor_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/donor_pipeline_${RUN_NAME}.err" "$donor_pipeline" "$donor_output_dir" "$new_param_file"
    echo "Finished running donor pipeline."
fi

# excision pipeline
if contains "excision" "${EDIT_STRATS[@]}"; then
    echo "Running excision pipeline..."
    excision_pipeline="$script_dir/edit_strategy_pipelines/excision_pipeline.sh"
    excision_output_dir=$output_dir"/excision"
    qsub -cwd -l mem_free=1G -l h_rt=08:00:00 -o "$project_root/logs/out/excision_pipeline_${RUN_NAME}.out" -e "$project_root/logs/err/excision_pipeline_${RUN_NAME}.err" "$excision_pipeline" "$excision_output_dir" "$new_param_file"
    echo "Finished running excision pipeline."
fi
