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

### Set Up Conda Environment
Create the conda environment from the repo root using the included `environment.yml` file:

```bash
conda env create -f environment.yml
```

This environment is named `excavate`. Activate it with:

```bash
conda activate excavate
```


### Data
We use data from several clinical and genomic databases including GTEx, ClinGen, and 1000 Genomes. Run the `data/data_dowloads.sh` script to ensure the necessary files are downloaded. Data storage requirements are listed in the `data` folder README.
```bash
bash ./data/data_downloads.sh
```

### Gene sets and running the pipeline
We constructructed the original D&D gene set using ClinGen Gene-Disease and Dosage Sensitivity summary information (`data/clingen`). `scripts/create_gene_set` contains the notebook used to select dominant and haplosufficient genes.

The main end-to-end pipeline entrypoint is:

```bash
bash ./scripts/run_edit_strategy_pipeline.sh
```

`run_edit_strategy_pipeline.sh` runs in two phases:

1. It runs the editing-strategy pipelines to generate per-strategy outputs.
2. After those finish, it launches `scripts/run_guide_analysis.sh` to perform the cross-strategy guide analysis.

Both phases use the same settings from `data/params/params.txt`.

You can also rerun just the guide-analysis phase on its own if the editing-strategy outputs already exist:

```bash
bash ./scripts/run_guide_analysis.sh <run_output_dir> <run_output_dir>/PARAMS/params.txt
```

## Running on an HPC
Code was originally intended to run on a high-performance compute cluster, and the current pipeline expects an SGE-style HPC environment. Included python scripts can additionally be run in a standalone format (outside of an HPC).

## Browser tracks
Code in `scripts/browser_tracks` is designed to generate viewable common variant tracks for genes of interest. Run
```bash
bash ./scripts/create_browser_tracks.sh
```
to generate bigBed files which can be uploaded to the UCSC Genome Browser TrackHub feature. Browser tracks for the current gene set can be found at the Github Repo (DnD_TrackHubs_Public)[https://github.com/gramey02/DnD_TrackHubs_Public] or at the session link https://genome.ucsc.edu/s/gramey02/MYH7_example_session.

## Outputs
Running the full pipeline produces information for each of the four gene editing strategies.
Alias names for these approaches in the results directory are:
- Exon disruption : "indels"
- Epigenetic silencing : "CRISPRoff"
- Splice site disruption : "donor/acceptor_base_edits"
- Excision : "excision"

In the results directory, example outputs include
- Per-gene editing strategy results in `results/<run_id>/<strategy>/`
- Number of common variants per gene and strategy in `results/<run_id>/<strategy>/ubiq_region_CommonVars/CommonVars_All_summary.txt`
- Common variant positions in `results/<run_id>/<strategy>/ubiq_region_CommonVars/CommonVars_All_dict.pkl`
- Summary tables (CSV) in `results/<run_id>/summary/`
- Browser track files (bigBed) in `results/browser_tracks/`

## Summary Outputs
You can generate summary information for your results by running the `scripts/make_summary_df/Creating_Master_DnD_DataFrame.ipynb`, followed by `scripts/make_summary_df/Final_df_formatting.ipynb`. These create the necessary outputs to pass into `scripts/figure_plotting/Figure_Plotting_Ntbk.ipynb`, which generates figures showing D&D gene properties and targetability.

## Gene Property Outputs
You can extract dominant pathogenic mutation counts and Human Phenotype Ontology (HPO) terms for each gene of interest. For mutation counts, run `scripts/clinvar_muts/ClinVar.ipynb`, and for ontology lists, run `scripts/hpo_terms/hpo_dominant_mapping.ipynb`.

## Citation
Please cite the coming preprint if you use this code.
