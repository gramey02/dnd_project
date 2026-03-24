#!/bin/bash
#$ -N find_excision_commonVars
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o logs/out/find_excision_commonVars.out
#$ -e logs/err/find_excision_commonVars.err

# Fail fast on errors, undefined variables, and pipeline failures.
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

output_dir="$1"
param_file="$2"
source "$param_file"
gene_info="$3"

script="$script_dir/find_excision_commonVars.py"
exon_file="$4"

python3 "$script" --output_dir "$output_dir" --gene_info "$gene_info" --af_limit "$AF_LIMIT" --af_file_dir "$AF_FILE_DIR" --exon_file "$exon_file" --nearby_gene_filter "$FILTER_OUT_NEARBY_GENES"
