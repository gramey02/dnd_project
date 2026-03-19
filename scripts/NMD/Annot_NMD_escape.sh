#!/bin/bash
#$ -N Annot_NMD_escape
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/Annot_NMD_escape.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/Annot_NMD_escape.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# parse input variables
# -----------------------------
output_dir="$1"
param_file="$2"
source "$param_file"
targetable_common_var_dict="$3"

# script to run NMD analysis
script="$script_dir/Annot_NMD_escape.py"
python3 "$script" --exon_file "$EXON_FILE_FOR_ANALYSIS" \
  --targetable_common_var_file "$targetable_common_var_dict" \
  --penultimate_rule "$PENULTIMATE_RULE" \
  --cds_rule "$CDS_RULE" \
  --exon_length_rule "$EXON_LENGTH_RULE" \
  --af_file_dir "$AF_FILE_DIR" \
  --output_dir "$output_dir"
