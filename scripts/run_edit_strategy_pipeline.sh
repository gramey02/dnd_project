#!/bin/bash
#$ -N run_edit_strategy_pipeline
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/run_edit_strategy_pipeline.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/run_edit_strategy_pipeline.err

# source params for the run
param_file="/wynton/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/params.txt"
source $param_file
echo $RUN_NAME

# ---- helper function: create a directory if missing, print a friendly message
mkdir_if_missing() {
  mkdir -p "$1"
}

# # make output directory for the run
output_dir=$OUTPUT_DIR$RUN_NAME
mkdir_if_missing $output_dir
param_dir=$output_dir"/PARAMS"
mkdir_if_missing $param_dir
cp $param_file $param_dir # save param file to the output directory for the run (so we can check later if needed)
new_param_file=$output_dir"/PARAMS/params.txt"

# run script to create the necessary directories for the results (right now, set up to have a directory per editing strategy)
make_dirs_script="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/make_dirs.sh"
bash $make_dirs_script $new_param_file

# run transcript filtering if param is set to true
if (( $FILTER_TRANSCRIPTS == 1 )); then
    mkdir -p "$output_dir/filtered_transcripts" # mkdir -p "$output_dir/filtered_transcripts/expression_filtered"
    # run transcript filtering
    filter_script="/wynton/protected/home/capra/gramey02/ConklinCollab/scripts/DN_CommonVars/filter_transcripts/filter_transcripts.sh" # include this in the params file eventually
    qsub -cwd -l mem_free=1G -l h_rt=01:00:00 -sync y $filter_script "$output_dir" "$new_param_file"
else
    # update
    printf 'EXON_FILE_FOR_ANALYSIS="%s"\n' "$ORIGINAL_EXON_FILE" >> "$new_param_file" # append KEY=VALUE
fi

# now call shell scripts to run sub-pipelines for each of the editing strategies and pass in the set parameters
# will parallelize these so they can run at the same time

# # indel pipeline
# echo "Running indel pipeline..."
# indel_pipeline="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/indel_pipeline.sh"
# indel_output_dir=$output_dir"/indels"
# qsub -cwd -l mem_free=1G -l h_rt=04:00:00 $indel_pipeline "$indel_output_dir" "$new_param_file"
# echo "Finished running indel pipeline."

# # crisproff pipeline
# echo "Running crisproff pipeline..."
# crisproff_pipeline="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/crisproff_pipeline.sh"
# crisproff_output_dir=$output_dir"/CRISPRoff"
# qsub -cwd -l mem_free=1G -l h_rt=04:00:00 $crisproff_pipeline "$crisproff_output_dir" "$new_param_file"
# echo "Finished running crisproff pipeline."

# # acceptor base edits pipeline
# echo "Running acceptor base edits pipeline..."
# acceptor_pipeline="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/acceptor_baseEdit_pipeline.sh"
# acceptor_output_dir=$output_dir"/acceptor_base_edits"
# qsub -cwd -l mem_free=1G -l h_rt=04:00:00 $acceptor_pipeline "$acceptor_output_dir" "$new_param_file"
# echo "Finished running acceptor pipeline."

# # donor base edits pipeline
# echo "Running donor base edits pipeline..."
# donor_pipeline="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/donor_baseEdit_pipeline.sh"
# donor_output_dir=$output_dir"/donor_base_edits"
# qsub -cwd -l mem_free=1G -l h_rt=04:00:00 $donor_pipeline "$donor_output_dir" "$new_param_file"
# echo "Finished running donor pipeline."

# excision pipeline
echo "Running excision pipeline..."
excision_pipeline="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/excision_pipeline.sh"
excision_output_dir=$output_dir"/excision"
qsub -cwd -l mem_free=1G -l h_rt=04:00:00 $excision_pipeline "$excision_output_dir" "$new_param_file"
echo "Finished running excision pipeline."


