# Dominant &amp; Dispensable (D&D) Gene Editing Project

## Overview
This code pertains to the D&D gene editing project, wherein we try to reduce the burden of getting bespoke gene editing therapies approved for D&D genes by therapeutically targeting common variation in the human genome.

## Quickstart
Clone repo using either
```bash
git clone git@github.com:gramey02/dnd_project.git
```
or
```bash
git clone https://github.com/gramey02/dnd_project.git
```

### 1 - Setup
__1a__. Start by setting up the conda environment.
Create the conda environment from the repo root using the included `environment.yml` file:

```bash
conda env create -f environment.yml
```

This environment is named `merged_env`. Activate it with:

```bash
conda activate merged_env
```
__1b__. Modify and run the setup script.
Navigate to the top level directory into which you cloned the repo, `dnd_project`. There is a setup script, `setup.sh`, that will replace remaining system-specific strings from where this code was originally run to strings specific to your system. Run this setup script prior to running analyses to ensure proper saving of files and email notifications about scripts running.
To run the setup script, navigate to the top level project directory. Run the following comands to print the working directory and make the setup script executable.
```bash
cd dnd_project
chmod +x setup.sh
pwd # should print something like '/Users/projects/dnd_project', or wherever you've cloned the repo to
```
Now, copy the working directory filepath that is printed. This will be the _first argument_ supplied to the setup script. The _second argument_ is your email. Note that you won't receive any emails unless you specify that as a parameter when running the analysis scripts later on.
Now you're ready to run the setup script! Run like so from the top-level project directory:
```bash
./setup.sh "<working_directory_name_you_just_printed>" "<your_preferred_email>"
```
Surrounding each argument in quotes is the safest option to avoid errors caused by spaces or other syntax issues.


### 2 - Data downloads
Now that we've modified the scripts to save to the right locations, we're ready to download data! We used data from several clinical and genomic databases including GTEx, ClinGen, and 1000 Genomes. Much of the data will be automatically downloaded upon cloning the repo, but files that are too large to store on GitHub will need to be downloaded from their web sources. Run the `data/data_dowloads.sh` script to ensure the necessary files are downloaded, which is set up with the appropriate download links. Please read the `data` folder's README to determine file size contraints for downloads and for more information on the versions of the data sources we used.
Run the data downloads script as follows, from a system or compute node that has internet access.
```bash
cd dnd_project
bash ./data/data_downloads.sh
```

### 3 - Creating gene sets
We constructructed the original D&D gene set using ClinGen Gene-Disease and Dosage Sensitivity summary information (`data/clingen`). `scripts/create_gene_set` contains the notebook used to select dominant and haplosufficient genes. We mapped ClinGen gene names to their HGNC approved names using the [HGNC biomart tool](https://biomart.genenames.org/).

In order to obtain information on transcripts, exon coordinates, and biotypes, we queried ENSEMBL. Gathering this information allows us to filter genes and transcripts by their biotype and determine where common genetic variants fall in the gene body. To query ENSEMBL for a new set of genes, one can run the script in `scripts/create_gene_set/get_ensembl_data.R` based on a list of gene names. Additionally, one can query directly from the [ENSEMBL Biomart website](https://useast.ensembl.org/info/data/biomart/index.html?).
Gene sets are saved under `data/dnd_hgnc` and their mapped ensembl data is saved under `data/dnd_ensembl_data.csv`.

### 4 - Running the pipeline
There are three main stages to the pipeline:
#### A. Finding genetic variants targetable by different gene editing strategies.
#### B. Finding gRNAs to target those variants.
#### C. Creating viewable UCSC Genome Browser tracks to visualize targetable variants.

A - To identify targetable genetic variants, the main entry point is the `run_edit_strategy_pipeline.sh` script. Run it using and SGE command like so:

```bash
cd dnd_project # IMPORTANT: navigate to top level project directory so relative filepaths don't break.
qsub -cwd -l mem_free=1G -l h_rt=00:10:00 ./scripts/run_edit_strategy_pipeline.sh
```
The `data/params/params.txt` file informs the parameters for this analysis. Please see the README in the params directory for information on what those parameters are and how to change them. The default parameters ensure the analysis is run for four different gene editing strategies, The output of this 'finding variants' analysis step The edit strategy pipeline outputs common variants in the 1000 Genomes population that are targetable by different CRISPR gene editing strategies, as well as the guide RNA (gRNA) sequences that would target each common variant.

Once the pipeline finishes running, users have the option to algorithmically rank the best gRNA sequences using a greedy algorithm approach.
The greedy algorithm prioritization can be run like so:
```bash
cd dnd_project # IMPORTANT: navigate to top level project directory so relative filepaths don't break.
qsub -cwd -l mem_free=1G -l h_rt=06:00:00 ./scripts/run_guide_analysis.sh
```
Prioritization of non-excision editing strategy guides will occur separately from excision editing strategy guide. Users have the option to specify in the params.txt file if they additionally want the gRNAs for each non-excision strategy to be considered separately or together during prioritization.

## Browser tracks
Code in `scripts/browser_tracks` is designed to generate viewable UCSC Genome Browser tracks for the targetable common variants of genes of interest. Run
```bash
cd dnd_project
qsub -cwd -l mem_free=1G -l h_rt=00:10:00 ./scripts/create_browser_tracks.sh
```
to generate bigBed files which can be uploaded to the UCSC Genome Browser TrackHub feature. Browser tracks for the current gene set can be found at the Github Repo [DnD_TrackHubs_Public](https://github.com/gramey02/DnD_TrackHubs_Public) or at the session link https://genome.ucsc.edu/s/gramey02/All_D%26Dgenes_w_filtering.

## Note on running scripts on an HPC
Code was originally intended to run on a high-performance compute (HPC) cluster, and the current pipeline expects an SGE-style HPC environment. Included python scripts can additionally be run in a standalone format (outside of an HPC).

## Outputs
Running the full pipeline produces information for each of the four gene editing strategies.
Alias names for the editing strategies in the results directory are:
- Exon disruption : "indels"
- Epigenetic silencing : "CRISPRoff"
- Splice site disruption : "donor_ or acceptor_base_edits"
- Excision : "excision"

In the results directory, example outputs include
- Per-gene editing strategy results in `results/<run_id>/<strategy>/`
- Number of common variants per gene and strategy in `results/<run_id>/<strategy>/ubiq_region_CommonVars/CommonVars_All_summary.txt`
- Common variant positions in `results/<run_id>/<strategy>/ubiq_region_CommonVars/CommonVars_All_dict.pkl`
- Summary tables (CSV) in `results/<run_id>/summary/`
- Browser track files (bigBed) in `results/browser_tracks/`

## Summary Outputs
You can generate summary information for your results by running the `scripts/make_summary_df/Creating_Master_DnD_DataFrame.ipynb` notebook, followed by `scripts/make_summary_df/Final_df_formatting.ipynb`. These create the necessary outputs to pass into `scripts/figure_plotting/Figure_Plotting_Ntbk.ipynb`, which generates figures showing D&D gene properties and targetability.

## Gene Property Outputs
You can extract dominant pathogenic mutation counts and Human Phenotype Ontology (HPO) terms for each gene of interest. For mutation counts, run `scripts/clinvar_muts/ClinVar.ipynb`, and for ontology lists, run `scripts/hpo_terms/hpo_dominant_mapping.ipynb`.

## Citation
Please cite the following [medRxiv preprint](https://www.medrxiv.org/content/10.64898/2026.03.26.26349431v1) when using code from this repository:
```bash
Ramey, G. D., Cowan, Q. T., Saxena, A. G., Macklin, B. L., Watry, H. L., Mei, S., ... & Capra, J. A. (2026). Leveraging human genetic variation to therapeutically target hundreds of genes with dominant & dispensable disease alleles. medRxiv, 2026-03.
```
