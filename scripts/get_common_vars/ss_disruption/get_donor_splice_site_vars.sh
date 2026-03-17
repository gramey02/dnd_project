#!/bin/bash
#$ -N get_donor_splice_site_vars
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/get_donor_splice_site_vars.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/get_donor_splice_site_vars.err

# parse input args
output_dir=$1
param_file=$2
source $param_file
exon_file=$3

script="/wynton/home/capra/gramey02/ConklinCollab/scripts/DN_BaseEditing/get_donor_splice_site_vars.py"

# Run the script below
python3 $script --exon_file $exon_file \
  --af_limit $AF_LIMIT \
  --editing_window_size $EDITING_WINDOW_SIZE \
  --donor_snp_region $DONOR_SNP_REGION \
  --output_dir $output_dir