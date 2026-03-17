#!/bin/bash
#$ -N find_excision_commonVars
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/find_excision_commonVars.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/find_excision_commonVars.err

output_dir=$1
param_file=$2
source $param_file
gene_info=$3

script="/wynton/home/capra/gramey02/ConklinCollab/scripts/pipeline_scripts/find_excision_commonVars.py"
exon_file=$4


python3 $script --output_dir $output_dir --gene_info $gene_info --af_limit $AF_LIMIT --exon_file $exon_file --nearby_gene_filt $FILTER_OUT_NEARBY_GENES

