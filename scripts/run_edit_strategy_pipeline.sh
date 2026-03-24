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

run_pipeline_async() {
  local strategy_name="$1"
  local pipeline_script="$2"
  local strategy_output_dir="$3"

  echo "Running ${strategy_name} pipeline..."
  bash "$pipeline_script" "$strategy_output_dir" "$new_param_file" &
  local pid=$!
  echo "Started ${strategy_name} pipeline with PID ${pid}."
  pipeline_pids+=("$pid")
  pipeline_names+=("$strategy_name")
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

exon_file="$EXON_FILE_FOR_ANALYSIS"
chrom_set_file="$output_dir/chromosomes/chrom_set.txt"
col=$(head -1 "$exon_file" | tr ',' '\n' | sed 's/^"//; s/"$//' | grep -nx "chromosome_name" | cut -d: -f1)
cut -d',' -f"$col" "$exon_file" | tail -n +2 | sed 's/^"//; s/"$//' | sort -u > "$chrom_set_file"

# now call shell scripts to run sub-pipelines for each of the requested editing strategies
pipeline_pids=()
pipeline_names=()

for strategy in "${EDIT_STRATS[@]}"; do
    case "$strategy" in
        indels)
            run_pipeline_async "indel" "$script_dir/edit_strategy_pipelines/indel_pipeline.sh" "$output_dir/indels"
            ;;
        CRISPRoff)
            run_pipeline_async "CRISPRoff" "$script_dir/edit_strategy_pipelines/crisproff_pipeline.sh" "$output_dir/CRISPRoff"
            ;;
        acceptor_base_edits)
            run_pipeline_async "acceptor" "$script_dir/edit_strategy_pipelines/acceptor_baseEdit_pipeline.sh" "$output_dir/acceptor_base_edits"
            ;;
        donor_base_edits)
            run_pipeline_async "donor" "$script_dir/edit_strategy_pipelines/donor_baseEdit_pipeline.sh" "$output_dir/donor_base_edits"
            ;;
        excision)
            run_pipeline_async "excision" "$script_dir/edit_strategy_pipelines/excision_pipeline.sh" "$output_dir/excision"
            ;;
        *)
            echo "Error: unsupported strategy '$strategy' in EDIT_STRATS." >&2
            exit 1
            ;;
    esac
done

pipeline_failures=0
for i in "${!pipeline_pids[@]}"; do
    pid="${pipeline_pids[$i]}"
    strategy_name="${pipeline_names[$i]}"
    if wait "$pid"; then
        echo "Finished running ${strategy_name} pipeline."
    else
        echo "Error: ${strategy_name} pipeline failed." >&2
        pipeline_failures=1
    fi
done

if (( pipeline_failures != 0 )); then
    exit 1
fi

# run cross-strategy guide analysis after all editing-strategy pipelines finish
if (( RUN_GUIDE_ANALYSIS == 1 )); then
    guide_analysis_script="$script_dir/run_guide_analysis.sh"
    echo "Running guide analysis..."
    bash "$guide_analysis_script" "$output_dir" "$new_param_file"
    echo "Finished running guide analysis."
else
    echo "Skipping guide analysis because RUN_GUIDE_ANALYSIS=$RUN_GUIDE_ANALYSIS."
fi
