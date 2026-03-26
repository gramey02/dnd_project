#!/bin/bash
#$ -N create_browser_tracks
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o logs/out/create_browser_tracks.out
#$ -e logs/err/create_browser_tracks.err

# -------------------------------
# Resolve project paths
# -------------------------------
# script_dir = where this script lives
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# project_root = repo root
project_root="$(cd "$script_dir/.." && pwd)"

# ensure log directories exist
mkdir -p "$project_root/logs/out" "$project_root/logs/err"

# -------------------------------
# Run promoter BED creation
# -------------------------------
bed_script="$script_dir/browser_tracks/create_gene_beds.sh"
targetable_genes_file=
num_targetable_genes=

qsub -cwd -sync y \
     -l mem_free=5G -l h_rt=02:00:00 \
     -t 1-$num_targetable_genes \
     -o "$project_root/logs/out/bed_creation.out" \
     -e "$project_root/logs/err/bed_creation.err" \
     "$bed_script"

echo "Finished promoter BED creation."

# -------------------------------
# Run bigBed conversion
# -------------------------------
bigbed_script="$script_dir/browser_tracks/create_bigbeds.sh"

qsub -cwd \
     -l mem_free=5G -l h_rt=02:00:00 \
     -o "$project_root/logs/out/bigbed_conversion.out" \
     -e "$project_root/logs/err/bigbed_conversion.err" \
     "$bigbed_script"

echo "Finished bigBed conversion."