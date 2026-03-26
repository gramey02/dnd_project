#!/bin/bash
#$ -N create_gene_beds
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

source "../../data/params/params.txt"
results_dir=$OUTPUT_DIR
run_name=$RUN_NAME
script="create_gene_beds.py"
gene_file="$results_dir/$run_name/summary_files/targetable_genes.txt"
ref_genome_fasta=$REF_GENOME_FASTA
bt_dir="../../data/browser_tracks/"
mkdir -p $bt_dir

# isolate one gene at a time for an array job
gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $gene_file)

# run script
python3 $script --results_dir $results_dir --run_name $run_name --gene $gene --ref_genome_fasta $ref_genome_fasta --bt_dir $bt_dir




