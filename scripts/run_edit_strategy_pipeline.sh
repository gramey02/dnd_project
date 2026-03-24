#!/bin/bash

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

# Resolve paths relative to this script rather than a fixed project location.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
execution_utils="$project_root/scripts/utils/execution_mode.sh"

# source params for the run
param_file="$project_root/data/params/params.txt"
source "$param_file"
source "$execution_utils"
echo "$RUN_NAME"
run_mode="$(get_run_mode)"
echo "Execution mode: $run_mode"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

# ---- helper function: create a directory if missing, print a friendly message
mkdir_if_missing() {
  mkdir -p "$1"
}

# Run a pipeline script locally with the output directory and params file.
run_pipeline() {
  local pipeline_script="$1"
  local strategy_output_dir="$2"
  bash "$pipeline_script" "$strategy_output_dir" "$new_param_file"
}

# # make output directory for the run
output_dir="${resolved_output_base}${RUN_NAME}"
mkdir_if_missing "$output_dir"
param_dir="$output_dir/PARAMS"
mkdir_if_missing "$param_dir"
cp "$param_file" "$param_dir" # save param file to the output directory for the run (so we can check later if needed)
new_param_file="$output_dir/PARAMS/params.txt"

# run script to create the necessary directories for the results (right now, set up to have a directory per editing strategy)
make_dirs_script="$script_dir/structure_result_dir/make_dirs.sh"
bash "$make_dirs_script" "$new_param_file"

# run transcript filtering if param is set to true
if (( $FILTER_TRANSCRIPTS == 1 )); then
    mkdir -p "$output_dir/filtered_transcripts" # mkdir -p "$output_dir/filtered_transcripts/expression_filtered"
    # run transcript filtering
    filter_script="$script_dir/filter_transcripts/filter_transcripts.sh" # include this in the params file eventually
    bash "$filter_script" "$output_dir" "$new_param_file"
else
    # update
    printf 'EXON_FILE_FOR_ANALYSIS="%s"\n' "$ORIGINAL_EXON_FILE" >> "$new_param_file" # append KEY=VALUE
fi

# now call shell scripts to run sub-pipelines for each of the editing strategies and pass in the set parameters
# will parallelize these so they can run at the same time

# indel pipeline
echo "Running indel pipeline..."
indel_pipeline="$script_dir/edit_strategy_pipelines/indel_pipeline.sh"
indel_output_dir="$output_dir/indels"
run_pipeline "$indel_pipeline" "$indel_output_dir"
echo "Finished running indel pipeline."

# crisproff pipeline
echo "Running crisproff pipeline..."
crisproff_pipeline="$script_dir/edit_strategy_pipelines/crisproff_pipeline.sh"
crisproff_output_dir="$output_dir/CRISPRoff"
run_pipeline "$crisproff_pipeline" "$crisproff_output_dir"
echo "Finished running crisproff pipeline."

# acceptor base edits pipeline
echo "Running acceptor base edits pipeline..."
acceptor_pipeline="$script_dir/edit_strategy_pipelines/acceptor_baseEdit_pipeline.sh"
acceptor_output_dir="$output_dir/acceptor_base_edits"
run_pipeline "$acceptor_pipeline" "$acceptor_output_dir"
echo "Finished running acceptor pipeline."

# donor base edits pipeline
echo "Running donor base edits pipeline..."
donor_pipeline="$script_dir/edit_strategy_pipelines/donor_baseEdit_pipeline.sh"
donor_output_dir="$output_dir/donor_base_edits"
run_pipeline "$donor_pipeline" "$donor_output_dir"
echo "Finished running donor pipeline."

# excision pipeline
echo "Running excision pipeline..."
excision_pipeline="$script_dir/edit_strategy_pipelines/excision_pipeline.sh"
excision_output_dir="$output_dir/excision"
run_pipeline "$excision_pipeline" "$excision_output_dir"
echo "Finished running excision pipeline."

# run cross-strategy guide analysis after all editing-strategy pipelines finish
guide_analysis_script="$script_dir/run_guide_analysis.sh"
echo "Running guide analysis..."
bash "$guide_analysis_script" "$output_dir" "$new_param_file"
echo "Finished running guide analysis."
