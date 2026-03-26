#!/bin/bash
#$ -N run_guide_analysis
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o logs/out/run_guide_analysis.out
#$ -e logs/err/run_guide_analysis.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/.." && pwd)"

output_dir="$1"
param_file="$2"

source "$param_file"
source "$execution_utils"

non_excision_guides_script="$script_dir/get_guides/non_excision_guides.sh"
echo "Running non-excision guide analysis..."
qsub -cwd -l mem_free=5G -l h_rt=02:00:00 "$non_excision_guides_script" "$output_dir" "$param_file"
echo "Finished non-excision guide analysis."

excision_guides_script="$script_dir/get_guides/excision_guides.sh"
echo "Running excision guide analysis..."
qsub -cwd -l mem_free=5G -l h_rt=02:00:00  "$excision_guides_script" "$output_dir" "$param_file"
echo "Finished excision guide analysis."
