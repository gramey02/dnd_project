# code to download data for dnd computations

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
data_dir="$(cd "$script_dir/.." && pwd)"

# downloading ClinGen Data
clingen_path="$data_dir/clingen"
mkdir -p "$clingen_path"
wget -O "$clingen_path/Clingen-Gene-Disease-Summary.csv" "https://search.clinicalgenome.org/kb/gene-validity/download" # gene-disease validity table
wget -O "$clingen_path/Clingen-Dosage-Sensitivity.csv" "https://search.clinicalgenome.org/kb/gene-dosage/download" # dosage sensitivity table

# downloading GENCODE data
gencode_path="$data_dir/GENCODE"
mkdir -p "$gencode_path"
wget -O "$gencode_path/gencode.v49.annotation.gtf.gz" https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz
gzip -dk "$gencode_path/gencode.v49.annotation.gtf.gz"

# downloading and unzipping GTEx v10 data
gtex_path="$data_dir/GTEx"
mkdir -p "$gtex_path"
wget -O "$gtex_path/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt" https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt # sample metadata
wget -O "$gtex_path/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz" "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz" # median TPMs per gene
wget -O "$gtex_path/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz" "https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz"
gzip -dkf "$gtex_path/GTEx_Analysis_v10_RNASeQCv2.4.2_gene_median_tpm.gct.gz"
gzip -dkf "$gtex_path/GTEx_Analysis_v10_RSEMv1.3.3_transcripts_tpm.txt.gz"

# download CpG Island data
cpg_path="$data_dir/cpg_islands"
mkdir -p "$cpg_path"
wget -O "$cpg_path/CpG_islands.txt.gz" "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cpgIslandExt.txt.gz"
if gzip -dc "$cpg_path/CpG_islands.txt.gz" | head -n 1 | grep -q "^#\|^chrom"; then
    gzip -dc "$cpg_path/CpG_islands.txt.gz" | tail -n +2 > "$cpg_path/CpG_islands_noHeader.txt"
else
    gzip -dc "$cpg_path/CpG_islands.txt.gz" > "$cpg_path/CpG_islands_noHeader.txt"
fi

# downloading 1000 Genomes data
onekg_dir="$data_dir/1KG"
mkdir -p "$onekg_dir"
wget -r -np -nH --cut-dirs=3 -R "index.html*" -P "$onekg_dir" https://hgdownload.soe.ucsc.edu/gbdb/hg38/1000Genomes/


# download hg38 fasta information
