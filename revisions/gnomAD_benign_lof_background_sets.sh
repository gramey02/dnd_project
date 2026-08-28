#!/bin/bash
#$ -N lof_background_recessive
#$ -t 1-1055 # for recessive there are 1055 genes # for dominant 262
#$ -cwd
#$ -o /wynton/home/capra/gramey02/dnd_project/revisions/logs/task_recessive_$TASK_ID.out
#$ -e /wynton/home/capra/gramey02/dnd_project/revisions/logs/task_recessive_$TASK_ID.err

# set up environment
module load CBI miniforge3 openjdk/11 aws-cli
conda activate dnscripts

dnd_file='/wynton/home/capra/gramey02/dnd_project/revisions/benign_lof_background_genes/recessive_genes_ensembl_mapped.csv'
ht_path='/wynton/home/capra/gramey02/dnd_project/revisions/hail_tables/gnomad_recessive_minus_dnd_pruned.ht'
output_dir='/wynton/home/capra/gramey02/dnd_project/revisions/results/recessive_minus_dnd/per_gene'

mkdir -p $output_dir
mkdir -p /wynton/home/capra/gramey02/dnd_project/revisions/checkpoints
mkdir -p /wynton/home/capra/gramey02/dnd_project/revisions/logs

# SGE_TASK_ID is 1-indexed; gene_spans.iloc is 0-indexed
gene_index=$((SGE_TASK_ID - 1))

script="/wynton/home/capra/gramey02/dnd_project/revisions/gnomAD_benign_lof_enrichment.py"
python3 "$script" --dnd_file "$dnd_file" --ht_path "$ht_path" --output_dir "$output_dir" --gene_index "$gene_index"