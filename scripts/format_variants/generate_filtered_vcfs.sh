#!/bin/bash
#$ -N generate_filtered_vcfs
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/ConklinCollab/scripts/out/generate_filtered_vcfs.out
#$ -e /wynton/home/capra/gramey02/ConklinCollab/scripts/err/generate_filtered_vcfs.err

module load CBI bcftools

# file that contains metadata on genes and coordinates
output_dir=$1
param_file=$2
source $param_file
gene_info=$3

# set up to run as an array job
cur_gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $gene_info)
cur_chrom=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' $gene_info)
# get current gene's common variant position file
common_vars_pos=$output_dir"/excavate/CommonVar_locs/${cur_gene}_CommonVar_locs.txt"

# output directory for the current analysis
#OUTPUT=
# filter the vcf on the current chromosome accordingly to these variants
output_vcf=$output_dir"/excavate/input_vcfs/${cur_gene}_CommonVar_filtered.vcf.gz"
# assign biallelic snp file based on current chromosome
biallelic_snps="$BIALLELIC_SNPS_DIR/TGP_chr${cur_chrom}_biallelicSNPs.vcf.gz"

# run filtering
bcftools view -T $common_vars_pos $biallelic_snps -Oz -o $output_vcf
# create index file
bcftools index -t $output_vcf



