#!/usr/bin/env bash
#$ -S /bin/bash
#$ -cwd

set -euo pipefail  # fail on errors, undefined vars, pipe issues

echo "Starting data downloads..."

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data_dir="$(cd "$script_dir/.." && pwd)"

########################################
# ClinGen
########################################
clingen_path="$data_dir/clingen"
mkdir -p "$clingen_path"

wget -nc -O "$clingen_path/Clingen-Gene-Disease-Summary.csv" \
  "https://search.clinicalgenome.org/kb/gene-validity/download"

wget -nc -O "$clingen_path/Clingen-Dosage-Sensitivity.csv" \
  "https://search.clinicalgenome.org/kb/gene-dosage/download"

########################################
# GENCODE
########################################
gencode_path="$data_dir/GENCODE"
mkdir -p "$gencode_path"

wget -nc -O "$gencode_path/gencode.v49.annotation.gtf.gz" \
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz

gzip -dkf "$gencode_path/gencode.v49.annotation.gtf.gz"

########################################
# GTEx
########################################
gtex_path="$data_dir/GTEx"
mkdir -p "$gtex_path"

wget -nc -O "$gtex_path/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt" \
  https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt

wget -nc -O "$gtex_path/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz" \
  https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz

wget -nc -O "$gtex_path/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz" \
  https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz

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
# hg38 reference (FIXED)
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
wget -nc -O "$ref_dir/hg38.chrom.sizes" \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

########################################

echo "All downloads complete!"