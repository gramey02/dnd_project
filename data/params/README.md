# Parameter file readme

This contains descriptions for each of the modifiable parameters listed in params.txt which are essential for outputs of the pipeline run.

## Transcript filtering parameters
* FILTER_TRANSCRIPTS - filter the input ensembl data based on transcript expression levels (accepted values = 0 or 1)
* GENE_EXPRESSION_PROP - TPM threshold at which we'll consider a gene to be 'expressed' in a tissue (float values accepted)
* EXPRESSION_PROP - Mean expression proportion of a transcript to be considered for downstream analyses (integer values accepted)
* KEEP_ALL_TRANSCRIPTS - if no transcripts for a gene pass the above thresholds, keep all transcripts for that gene. If False, keeps only the canonical transcript for the gene.
* TRANSCRIPT_TPM_FILE - file with per transcript TPMs per sample per tissue
* SAMPLE_ATTRIBUTES_FILE - file that maps tissue name to sample ID
* GENE_MEDIAN_TPMS_FILE - file with median transcript TPM per tissue across samples
* COLNAMES_FILE - file that maps transcript names easily
* TISSUE_MAP_FILE - file that maps tissue names with underscores to tissue names with spaces
* ORIGNINAL_EXON_FILE - file with exon information and transcript biotype information per gene of interest, from ENSEMBL

## Common variants identification parameters
* AF_LIMIT - allele frequency threshold that we will consider 'common'. Values >= to this value are considered passing.
* AF_FILE_DIR - directory containing 1000Genomes variant allele frequencies
* BIALLELIC_SNPS_DIR - directory containing full vcfs of biallelic snps in 1000Genomes
* PHASED_1000G_VCF_DIR= - directory containing all original phased vcfs from 1000Genomes
* USE_UBIQUITOUS_REGIONS - _depracated and unused_
* USE_ALL_REGIONS - _depracated and unused_
* REQUIRED_TRANSCRIPT_PROPORTION - _depracated and unused_

## Nonsense mediated decay analysis parameters
* RUN_NMD - run nonsense mediated decay analysis, or if all variants should be passed to next step (accepts 0 or 1)
* PENULTIMATE_RULE - number of base pairs to consider for penultimate rule. See aenmd publication for number of bp in penultimate exon that escape NMD
* CDS_RULE - number of base pairs to consider for coding sequence rule. See aenmd publication for coding sequence rule on which variants escape NMD
* EXON_LENGTH_RULE - number of base pairs to consider for exon length rule. See aenmd publication for exon lenght rule on which variants escape NMD

## Epigenetic silencing analysis parameters
* PROMOTER_UD - number of base pairs upstream of transcription start site (TSS) to consider promoter start
* PROMOTER_DD - number of base pairs upstream of TSS to consider promoter end
* INTERSECTION - minimum acceptable intersection size to consider a genomic region overlapping with a CpG Island
* USE_ISLANDS - whether or not CpG Island overlap should be used to filter promoter regions, or if just GC content should be used (see next param)
* GC_THRESH - minimum percentage threshold of GC content that a promoter should meet to be considered a CRISPR-mediated epigenetically silencable target
* CPG_FILE - file with CpG island genomic intervals

## Base editing strategy parameters
Note that we have separated donor splice sites and acceptor splice sites, but later combine their results.
* EDITING_WINDOW_SIZE - optimal editing window size of base editors (e.g., ABEs, CBEs)
* DONOR_SNP_REGION - Window in which we should look for common genetic variants nearby splice donor sites
* ACCEPTOR_SNP_REGION - Window in which we should look for common genetic variants nearby splice acceptor sites

## EXCAVATE parameters
* REF_GENOME_FASTA - filepath to hg38 human reference genome fasta file
* CHROM_FASTA_FILEPATH - filepath to human reference genome fastas separated by chromosome
* EXCAVATE_SCRIPT - filepath to EXCAVATE entrypoint script

## Excision pipeline parameters
* UPSTREAM_EXCISION_DIST - base pair distance upstream of a gene to search for excision variants
* DOWNSTREAM_EXCISION_DIST - base pair distance downstream of a gene to search for excision variants
* GENE_BODY_FILE - filepath to GENCODE file with gene body coordinates
* EXCISE_ENTIRE_GENE - if true, pipeline considers only genetic variants that flank the entire gene for excision, rather than ones that flank indivdual exons
* NUM_SAMPLES=2548 - number of samples in a population. Commonly the 1000Genomes population
* FILTER_OUT_NEARBY_GENES - if true, pipeline restricts excision window by removing nearby protein coding gene coordinates from excision window

## All edit strategy parameters
* EDIT_STRATS - tells the pipeline wrapper which edit strategy sub-pipelines to run. Default is exon disruption, splice site disruption, epigenetic silencing, and excision

## Parameters for storing project outputs
* PROJECT_ROOT - the location of dnd_project folder (e.g., '/User/projects/dnd_project')
* RUN_NAME - name of the pipeline run. Defaults to current date
* OUTPUT_DIR - location of results directory
* SUMMARY_FILE - create a summary file distilling the information obtained for each editing strategy (e.g., how many genes are captured by each strategy, how many heterozygous individuals are captured, etc.)
