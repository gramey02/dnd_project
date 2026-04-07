#!/bin/bash
#$ -N run_excavate
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# activate excavate environment
module load CBI miniforge3 bcftools
conda activate excavate

# load input vars
output_dir=$1
param_file=$2
source $param_file
gene_info=$3
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# set variables for inputs and outputs
ref_genome_fasta="$project_root/$REF_GENOME_FASTA"
gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $gene_info)
chrom=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' $gene_info)
gene_coords=$(awk -v row=$SGE_TASK_ID 'NR == row {print $5}' $gene_info)
common_var_vcf=$output_dir"/excavate/input_vcfs/${gene}_CommonVar_filtered.vcf.gz"
base="${CHROM_FASTA_FILEPATH%/}"
chrom_fasta="$project_root/$base/hg38_${chrom}.fa"

# troubleshooting
# echo $ref_genom_fasta
# echo $gene
# echo $chrom
# echo $gene_coords
# echo $common_var_vcf
# echo $chrom_fasta

# make output file directory below
output_files=$output_dir"/excavate/excavate_outputs/${gene}_output"
if [ ! -d "$output_files" ]; then
    mkdir -p "$output_files"
fi

echo $output_files

# # run excavate for each gene below
# python3 "$project_root/$EXCAVATE_SCRIPT" generate $common_var_vcf population $chrom_fasta $ref_genome_fasta $gene_coords --cas SpCas9 --summary -o $output_files
