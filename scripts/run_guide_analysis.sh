#!/bin/bash
#$ -N run_guide_analysis
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o logs/out/run_guide_analysis.out
#$ -e logs/err/run_guide_analysis.err

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"
execution_utils="$project_root/scripts/utils/execution_mode.sh"

output_dir="$1"
param_file="$2"

source "$param_file"
source "$execution_utils"

run_mode="$(get_run_mode)"
echo "Running guide analysis in $run_mode mode..."

non_excision_guides_script="$script_dir/get_guides/non_excision_guides.sh"
echo "Running non-excision guide analysis..."
bash "$non_excision_guides_script" "$output_dir" "$param_file"
echo "Finished non-excision guide analysis."

excision_guides_script="$script_dir/get_guides/excision_guides.sh"
echo "Running excision guide analysis..."
bash "$excision_guides_script" "$output_dir" "$param_file"
echo "Finished excision guide analysis."
