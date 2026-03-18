# code to download data for dnd computations

# downloading ClinGen Data
wget https://search.clinicalgenome.org/kb/gene-validity/download # gene-disease validity table
wget https://search.clinicalgenome.org/kb/gene-dosage/download # dosage sensitivity table

# downloading GENCODE data
gencode_path=".././GENCODE/"
mkdir -p $gencode_path
wget -O "$gencode_path/gencode.v49.annotation.gtf.gz" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
gzip -dk "$gencode_path/gencode.v49.annotation.gtf.gz"

# downloading and unzipping GTEX v10 data
gtex_path=".././GTEx_v10/"
mkdir -p $gtex_path
wget -O "$gtex_path/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt" https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt # sample metadata
wget -O "$gtex_path/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz" "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz" # median TPMs per gene
wget -O "$gtex_path/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz" "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz"
gzip -dk "$gtex_path/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
gzip -dk "$gtex_path/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz"

# downloading 1000 Genomes data


# download hg38 information

