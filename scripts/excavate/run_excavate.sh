#!/bin/bash
#$ -N run_excavate
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../.." && pwd)"

# activate excavate environment
module load CBI miniforge3 bcftools
conda activate excavate

# load input vars
output_dir=$1
param_file=$2
source $param_file
gene_info=$3


# load input vars
output_dir=$1
param_file=$2
source $param_file
gene_info=$3

# set variables for inputs and outputs
ref_genome_fasta=$REF_GENOME_FASTA
gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $gene_info)
chrom=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' $gene_info)
gene_coords=$(awk -v row=$SGE_TASK_ID 'NR == row {print $5}' $gene_info)
common_var_vcf=$output_dir"/excavate/input_vcfs/${gene}_CommonVar_filtered.vcf.gz"
chrom_fasta=$CHROM_FASTA_FILEPATH"hg38_${chrom}.fa"

# make output file directory below
output_files=$output_dir"/excavate/excavate_outputs/${gene}_output"
if [ ! -d "$output_files" ]; then
    mkdir -p "$output_files"
    echo "Directory '$output_files' created."
else
    echo "Directory '$output_files' already exists."
fi

# run excavate for each gene below
python3 $EXCAVATE_SCRIPT generate $common_var_vcf population $chrom_fasta $ref_genome_fasta $gene_coords --cas SpCas9 --summary -o $output_files
