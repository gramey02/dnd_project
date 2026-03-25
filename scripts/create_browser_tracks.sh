#!/bin/bash
#$ -N create_browser_tracks
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ../logs/out/create_browser_tracks.out
#$ -e ../logs/err/create_browser_tracks.err

# run bed creation script
qsub -cwd -sync y -l mem_free=5G -l h_rt=02:00:00 ./browser_tracks/create_gene_beds.sh

# run bigBed conversion script
qsub -cwd -l mem_free=5G -l h_rt=02:00:00 ./browser_tracks/create_bigbeds.sh
