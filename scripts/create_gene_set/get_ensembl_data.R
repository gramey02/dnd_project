# install packages if missing, then load libraries
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("biomaRt", quietly = TRUE)) {
  BiocManager::install("biomaRt", ask = FALSE, update = FALSE)
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}

# load packages
library(biomaRt)
library(dplyr)

# resolve paths relative to this script
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- dirname(normalizePath(script_path))
project_root <- normalizePath(file.path(script_dir, "..", ".."))

# load ensembl biomaRt
ensembl<-useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")

# load data that contains genes of interest
input_file <- file.path(project_root, "data", "dnd_hgnc", "dhs_hgnc_mapped.csv")
df<-read.csv(input_file)

# queries
# need to split them because some attributes lie on different attribute pages
res1<-getBM(attributes=c('hgnc_symbol','ensembl_gene_id','ensembl_transcript_id','chromosome_name','start_position','end_position','strand', 'transcript_start', 'transcript_end',
'transcription_start_site','transcript_length','transcript_is_canonical','transcript_tsl','gene_biotype','transcript_biotype','ensembl_peptide_id'),
filters = 'ensembl_gene_id',values=df$Ensembl.gene.ID,mart=ensembl)

res2<-getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','exon_chrom_start','exon_chrom_end','ensembl_exon_id','is_constitutive','rank','genomic_coding_start','genomic_coding_end','cds_start','cds_end'),
filters = 'ensembl_gene_id',values=df$Ensembl.gene.ID,mart=ensembl)

# merge and save the data
merged<-merge(res1,res2,by=c('ensembl_gene_id','ensembl_transcript_id'),all.x = TRUE, all.y=TRUE)
output_file <- file.path(project_root, "data", "dnd_ensembl", "dnd_ensembl_data.csv")
write.csv(merged, output_file, row.names = FALSE)


# ==================================================================================================
# in case you're running it straight from the multi-symbol checker step, you can do the following steps:

# # load libraries
# library(biomaRt)
# library(dplyr)

# # load ensembl biomaRt
# ensembl<-useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
# df<-read.csv('/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_11_13/dhs_hgnc_mapped_2025-11-13.csv')
# df_filt <- df %>% filter(Match.type=='Approved symbol')

# res1<-getBM(attributes=c('hgnc_symbol','ensembl_gene_id','ensembl_transcript_id','chromosome_name','start_position','end_position',
# 'strand', 'transcript_start', 'transcript_end','transcription_start_site','transcript_length','transcript_is_canonical','transcript_tsl',
# 'gene_biotype','transcript_biotype','ensembl_peptide_id'),filters = 'hgnc_symbol',values=df_filt$Approved.symbol,mart=ensembl)

# # it's the additional chromosome annotations that add to the ensembl identifiers, so standardize these here
# res1_filt<-res1 %>% filter(chromosome_name %in% c('1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X'))

# res2<-getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id','exon_chrom_start','exon_chrom_end','ensembl_exon_id',
# 'is_constitutive','rank','genomic_coding_start','genomic_coding_end','cds_start','cds_end'),
# filters = 'ensembl_gene_id',values=res1_filt$ensembl_gene_id,mart=ensembl)

# # save
# merged<-merge(res1_filt,res2,by=c('ensembl_gene_id','ensembl_transcript_id'),all.x = TRUE, all.y=TRUE)
# write.csv(merged,'/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_11_13/dHS_exonInfo_2025_11_13.csv')
