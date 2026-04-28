#!/bin/bash
#$ -N create_gene_beds
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input args
param_file="$1"
gene_file="$2"
output_dir="$3"

# load variables from input args
source $param_file
project_root="$PROJECT_ROOT"
results_dir="$PROJECT_ROOT/$OUTPUT_DIR"
script_dir="$project_root/scripts"
run_name=$RUN_NAME

# set scripts
script="$script_dir/browser_tracks/create_gene_beds.py"
bt_dir="$output_dir/summary_files/browser_tracks"
exon_file=$ORIGINAL_EXON_FILE
mkdir -p $bt_dir

# isolate one gene at a time for an array job
gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $gene_file)
# get current chromosome for the gene
cur_chrom=$(awk -F',' -v gene="$gene" 'NR==1{for(i=1;i<=NF;i++){gsub(/"/,"",$i); if($i=="hgnc_symbol")g=i; if($i=="chromosome_name")c=i}} NR>1 && $g==gene{val=$c; gsub(/"/,"",val); print val; exit}' "$PROJECT_ROOT/$ORIGINAL_EXON_FILE")
rsID_fp="$project_root/$RSID_MAPS/gnomAD_v4.1.0_PASS_SNVs_chr${cur_chrom}_pos_rsIDs.pkl.gz"


# run script
python3 $script --results_dir $results_dir --run_name $run_name --gene $gene \
        --ref_genome_fasta $REF_GENOME_FASTA --bt_dir $bt_dir --exon_file $exon_file \
        --sample_map $SAMPLE_MAP --rsID_fp $rsID_fp

