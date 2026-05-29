#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

echo "Starting data downloads..."

param_file="/wynton/home/capra/gramey02/dnd_project/data/params/params.txt"
source "$param_file"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"
data_dir="$project_root/data"

# activate conda environment, since we'll need bcftools
conda activate merged_env

########################################
# ClinGen
########################################
# # uncomment to download from the source
# # for now, it's already included in the original repo
# clingen_path="$data_dir/clingen"
# mkdir -p "$clingen_path"

# wget -nc -O "$clingen_path/Clingen-Gene-Disease-Summary.csv" \
#   "https://search.clinicalgenome.org/kb/gene-validity/download"

# wget -nc -O "$clingen_path/Clingen-Dosage-Sensitivity.csv" \
#   "https://search.clinicalgenome.org/kb/gene-dosage/download"



########################################
# GENCODE
########################################
gencode_path="$data_dir/GENCODE"
mkdir -p "$gencode_path"

# # uncomment to download from the source
# # for now, it's already included in the original repo
# wget -nc -O "$gencode_path/gencode.v49.annotation.gtf.gz" \
#   https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz

# unzip file
gzip -dkf "$gencode_path/gencode.v49.annotation.gtf.gz"



########################################
# GTEx
########################################
gtex_path="$data_dir/GTEx"
mkdir -p "$gtex_path"

# # uncomment if downloading again from the source
# # for now, these are already included in the original repo
# wget -nc -O "$gtex_path/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt" \
#   https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt
# wget -nc -O "$gtex_path/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz" \
#   https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz

# this file is too large to store on Github, so needs to be downloaded
wget -nc -O "$gtex_path/dhsGenes_GTEx_transcript_TPMs.txt" \
  https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz
# unzip the file
gzip -dkf "$gtex_path/"*.gz



########################################
# CpG Islands
########################################
cpg_path="$data_dir/cpg_islands"
mkdir -p "$cpg_path"

wget -nc -O "$cpg_path/CpG_islands.txt.gz" \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz

if gzip -dc "$cpg_path/CpG_islands.txt.gz" | head -n 1 | grep -q "^#\|^chrom"; then
    gzip -dc "$cpg_path/CpG_islands.txt.gz" | tail -n +2 > "$cpg_path/CpG_islands_noHeader.txt"
else
    gzip -dc "$cpg_path/CpG_islands.txt.gz" > "$cpg_path/CpG_islands_noHeader.txt"
fi


########################################
# 1000 Genomes (UCSC mirror)
########################################
onekg_dir="$data_dir/1KG"
mkdir -p "$onekg_dir"

wget -r -np -nH --cut-dirs=3 -R "index.html*" \
  -P "$onekg_dir" \
  https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/



########################################
# hg38 reference
########################################
ref_dir="$data_dir/reference_genomes"
mkdir -p "$ref_dir"

# Option A (recommended): full genome fasta
wget -nc -O "$ref_dir/hg38.fa.gz" \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gzip -dkf "$ref_dir/hg38.fa.gz"

# Option B (if you truly want per-chromosome files, uncomment):
# wget -r -np -nH --cut-dirs=4 -R "index.html*" \
#   -P "$ref_dir" \
#   ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/



########################################
# Chrom sizes (FIXED PATH)
########################################
# uncomment if downloading from the source
# # for now, this file is already included in the repository
# wget -nc -O "$ref_dir/hg38.chrom.sizes" \
#   https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

########################################

# unzip files that need to be unzipped
gzip -dk "$project_root/data/biallelic_afs/TGP*"
gzip -dk "$project_root/data/GENCODE/gencode.v49.annotation.gtf.gz"



# make the biallelic vcf files
mkdir -p "$data_dir/biallelic_vcfs"

for f in "$data_dir/1KG/ALL.chr*.vcf.gz"; do
    chr=$(basename "$f" | cut -d'.' -f2 | sed 's/chr//')
    out="data/biallelic_vcfs/TGP_chr${chr}_biallelicSNPs.vcf.gz"
    bcftools view \
        -m2 -M2 \
        -v snps \
        -Oz \
        -o "$out" \
        "$f"
    bcftools index "$out"
done

echo "All downloads complete!"