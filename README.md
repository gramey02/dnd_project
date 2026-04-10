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
1a. Start by setting up the conda environment
Create the conda environment from the repo root using the included `environment.yml` file:

```bash
conda env create -f environment.yml
```

This environment is named `excavate`. Activate it with:

```bash
conda activate excavate
```
2b. Modify and run the setup script
Navigate to `scripts/setup.sh`.
Replace "<your_email_here>" with your preferred email to receive notifications about the pipeline run.
Replace "<dnd_project_directory>" with the full directory into which you cloned the dnd repo. You can find this full directory string by navigating to dnd_project and printing the working directory:
```bash
cd dnd_project
pwd # should print something like '/Users/projects/dnd_project', or wherever you've cloned the repo to
```
Next, run the setup script:
```bash
cd scripts
bash ./setup.sh
```
This script finds and replaces filepaths throughout scripts that will be specific to your computer, so results can be stored in the right location and so the right filepaths can be found while the pipeline is running.


### 2 - Data downloads
We use data from several clinical and genomic databases including GTEx, ClinGen, and 1000 Genomes. Run the `data/data_dowloads.sh` script to ensure the necessary files are downloaded. Data storage requirements are listed in the `data` folder README.
```bash
bash ./data/data_downloads.sh
```

### 3 - Creating gene sets and running the pipeline
We constructructed the original D&D gene set using ClinGen Gene-Disease and Dosage Sensitivity summary information (`data/clingen`). `scripts/create_gene_set` contains the notebook used to select dominant and haplosufficient genes.

The main end-to-end pipeline entrypoint is the run_edit_strategy_pipeline.sh script. Run it using and SGE command like so:

```bash
cd dnd_project # IMPORTANT: navigate to top level project directory so relative filepaths don't break.
qsub -cwd -l mem_free=1G -l h_rt=00:10:00 ./scripts/run_edit_strategy_pipeline.sh
```
The edit strategy pipeline outputs common variants in the 1000 Genomes population that are targetable by different CRISPR gene editing strategies, as well as the guide RNA (gRNA) sequences that would target each common variant.

Once the pipeline finishes running, users have the option to algorithmically rank the best gRNA sequences using a greedy algorithm approach.
The greedy algorithm prioritization can be run like so:
```bash
cd dnd_project # IMPORTANT: navigate to top level project directory so relative filepaths don't break.
qsub -cwd -l mem_free=1G -l h_rt=00:10:00 ./scripts/run_guide_analysis.sh
```
Prioritization of non-excision editing strategy guides will occur separately from excision editing strategy guide. Users have the option to specify in the params.txt file if they additionally want the gRNAs for each non-excision strategy to be considered separately or together during prioritization.

## Browser tracks
Code in `scripts/browser_tracks` is designed to generate viewable UCSC Genome Browser tracks for the targetable common variants of genes of interest. Run
```bash
cd dnd_project
qsub -cwd -l mem_free=1G -l h_rt=00:10:00 ./scripts/create_browser_tracks.sh
```
to generate bigBed files which can be uploaded to the UCSC Genome Browser TrackHub feature. Browser tracks for the current gene set can be found at the Github Repo [DnD_TrackHubs_Public](https://github.com/gramey02/DnD_TrackHubs_Public) or at the session link https://genome.ucsc.edu/s/gramey02/MYH7_example_session.

## Note on running scripts on an HPC
Code was originally intended to run on a high-performance compute cluster, and the current pipeline expects an SGE-style HPC environment. Included python scripts can additionally be run in a standalone format (outside of an HPC). Adaption of the code to run locally and on Slurm-based clusters is in development.

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
