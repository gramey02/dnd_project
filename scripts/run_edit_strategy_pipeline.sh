#!/bin/bash

set -euo pipefail

# Resolve paths relative to this script rather than a fixed project location.
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

# source params for the run
param_file="$project_root/data/params/params.txt"
source "$param_file"
echo "$RUN_NAME"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

# ---- helper function: create a directory if missing, print a friendly message
mkdir_if_missing() {
  mkdir -p "$1"
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

# Additional editing-strategy pipelines can be invoked here with plain
# `bash "$pipeline_script" "$strategy_output_dir" "$new_param_file"` calls.

# excision pipeline
echo "Running excision pipeline..."
excision_pipeline="$script_dir/edit_strategy_pipelines/excision_pipeline.sh"
excision_output_dir="$output_dir/excision"
bash "$excision_pipeline" "$excision_output_dir" "$new_param_file"
echo "Finished running excision pipeline."
