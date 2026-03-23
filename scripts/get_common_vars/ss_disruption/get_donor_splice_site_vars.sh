#!/bin/bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# parse input args
output_dir="$1"
param_file="$2"
source "$param_file"
exon_file="$3"

script="$script_dir/get_donor_splice_site_vars.py"

# Run the script below
python3 "$script" --exon_file "$exon_file" \
  --af_limit "$AF_LIMIT" \
  --af_file_dir "$AF_FILE_DIR" \
  --editing_window_size "$EDITING_WINDOW_SIZE" \
  --donor_snp_region "$DONOR_SNP_REGION" \
  --output_dir "$output_dir"
