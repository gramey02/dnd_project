#!/bin/bash
#$ -N promoter_common_vars
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/promoter_common_vars.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/promoter_common_vars.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# parse input params
output_dir="$1"
param_file="$2"
source "$param_file"
exon_file="$3"

script="$script_dir/promoter_common_vars.py"
python3 "$script" --cpg_file "$CPG_FILE" \
  --gt_file "$exon_file" \
  --intersection "$INTERSECTION" \
  --promoter_ud "$PROMOTER_UD" \
  --promoter_dd "$PROMOTER_DD" \
  --output_dir "$output_dir" \
  --af_limit "$AF_LIMIT" \
  --gc_threshold "$GC_THRESH" \
  --use_islands "$USE_ISLANDS" \
  --ref_genome_fasta "$REF_GENOME_FASTA"
