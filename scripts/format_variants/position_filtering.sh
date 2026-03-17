#!/bin/bash
#$ -N position_filtering
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ~/ConklinCollab/scripts/out/position_filtering.out
#$ -e ~/ConklinCollab/scripts/err/position_filtering.err

module load CBI bcftools

# load input arguments
output_dir=$1
param_file=$2
source $param_file
gene_info=$3


# separate out inputs and outputs
chrom=$(awk -v row=$SGE_TASK_ID 'NR == row {print $2}' $gene_info)
gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $gene_info)
positions=$output_dir"/excavate/Guide_locs/${gene}_Guide_locs.txt"
input_vcf="/wynton/group/capra/data/wynton_databases/1000_genomes/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
filtered_vcf=$output_dir"/excavate/Guide_filtered_vcfs/${gene}_guide_filtered.vcf.gz"

# filter with bcftools
bcftools view -T $positions $input_vcf -Oz -o $filtered_vcf
# unzip so you can run analysis on the file later
gzip -dk $filtered_vcf









# non-parallelized
# # test:
# gene_info="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/indel_PAM_filtered_vars/dhs_small.txt"
# # full:
# #gene_info="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/indel_PAM_filtered_vars/dhs_indel_PAM_2025_04_22_noIDX.txt"

# edit_strategy="indel"

# while read gene chrom start_pos end_pos gene_coords; do
#     echo "Started ${gene} filtering..."
#     # assign variables for the current gene
#     positions="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/${edit_strategy}_PAM_filtered_vars/${gene}_${edit_strategy}_PAM_positions.txt"
#     input_vcf="/wynton/group/capra/data/wynton_databases/1000_genomes/release/20190312_biallelic_SNV_and_INDEL/ALL.chr${chrom}.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz"
#     filtered_vcf="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/${edit_strategy}_PAM_filtered_vars/${gene}.vcf.gz"
#     bcftools view -T $positions $input_vcf -Oz -o $filtered_vcf
#     gzip -dk $filtered_vcf
#     echo "Finished ${gene} filtering."
# done < $gene_info
