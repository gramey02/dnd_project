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
Create the conda environment from the repo root using the included [environment.yml](/wynton/home/capra/gramey02/dnd_project/environment.yml) file:

```bash
conda env create -f environment.yml
```

This environment is named `excavate`. Activate it with:

```bash
conda activate excavate
```


### Data downloads
We use data from several clinical and genomic databases including GTEx, ClinGen, and 1000 Genomes. Run the data_dowloads.sh script to ensure the necessary files are downloaded.
```bash

```

### Gene sets and running the pipeline
We constructructed the original D&D gene set using ClinGen Gene-Disease and Dosage Sensitivity summary information (data/clingen). `scripts/create_gene_set` contains the notebook used to select dominant and haplosufficient genes, and ` 


## Running on an HPC
Code was originally intended to run on an HPC. Runtime for the entire pipeline excluding jupyter notebooks for data visualization ranged from 4-5 hours for the original ~600 D&D genes of interest. No GPUs were necessary, but given intensity of the greedy algorithm gRNA prioritization approach, they are recommended.

## Summary Outputs
You can generate summary information for your results by running the `scripts/make_summary_df/Creating_Master_DnD_DataFrame.py`, followed by `scripts/make_summary_df/Final_df_formatting.py`. These create the necessary outputs to pass into `scripts/figure_plotting/Figure_Plotting_Ntbk.ipynb`, which generates figures showing D&D gene properties and targetability.

## Citation
Please cite the coming preprint if you use this code.
