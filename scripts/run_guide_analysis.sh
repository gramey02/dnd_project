#!/bin/bash
#$ -N run_guide_analysis
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o logs/out/run_guide_analysis.out
#$ -e logs/err/run_guide_analysis.err

# parse input args
project_root="/wynton/home/capra/gramey02/dnd_project"
script_dir="$project_root/scripts"

# source params for the run
param_file="$project_root/data/params/params.txt"
source "$param_file"

# load parameters from sources file
output_dir="${OUTPUT_DIR}${RUN_NAME}"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# run the non-excision guides analysis, in which one guide is selected upon each iteration
non_excision_guides_script="$script_dir/get_guides/non_excision_guides.sh"
echo "Running non-excision guide analysis..."
qsub -cwd -l mem_free=5G -l h_rt=02:00:00 -o "$project_root/logs/out/non_excision_guides.out" -e "$project_root/logs/err/non_excision_guides.err" "$non_excision_guides_script" "$project_root/$output_dir" "$param_file"

# run the excision guides analysis, in which one or two guides is selected upon each iteration (max 50 iterations)
excision_guides_script="$script_dir/get_guides/excision_guides.sh"
echo "Running excision guide analysis..."
qsub -cwd -l mem_free=5G -l h_rt=02:00:00  -o "$project_root/logs/out/excision_guides.out" -e "$project_root/logs/err/excision_guides.err" "$excision_guides_script" "$project_root/$output_dir" "$param_file"
