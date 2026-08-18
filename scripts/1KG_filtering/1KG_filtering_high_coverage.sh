#!/bin/bash
#$ -S /bin/bash
#$ -N biallelic_afs
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/home/capra/gramey02/dnd_project/scripts/1KG_filtering/logs/biallelic_afs_chrX.out
#$ -e /wynton/home/capra/gramey02/dnd_project/scripts/1KG_filtering/logs/biallelic_afs_chrX.err
##$ -t 1-23
#$ -l mem_free=8G
#$ -l h_rt=04:00:00
#
# SGE array job: one task per chromosome (1-22, then X as task 23)
# For each chromosome:
#   1. Filter CCDG high-coverage 1000G VCF to biallelic SNPs only
#   2. Subset to the 2,504 unrelated samples (drops the 698 additional related samples)
#   3. Recompute AC/AN/AF (overall + per-superpopulation) on the subsetted genotypes,
#      using -S with a sample->superpopulation mapping file (fill-tags generates the
#      AC_<GROUP>/AN_<GROUP>/AF_<GROUP> tags automatically -- they can't be requested
#      directly by name via -t)
#   4. Rename chromosome from chrN -> N (both header ##contig lines and data lines)
#   5. Save filtered+subsetted+renamed VCF to biallelic_vcfs_high_coverage/
#   6. Extract CHROM,POS,REF,ALT,AC,AN,AF,AF_EAS,AF_EUR,AF_AFR,AF_AMR,AF_SAS (no header)
#   7. Save summary to biallelic_afs_high_coverage/
#
# IMPORTANT: run `mkdir -p logs` in the target log directory before qsub-ing, or
# SGE will fail to write the -o/-e log files and the job will error immediately.
#
# Usage: qsub 1KG_filtering_high_coverage.sh

set -euo pipefail  # exit on error, undefined var, or failed pipe -- catches silent failures early
module load CBI bcftools

# ---- EDIT THESE PATHS ----
INPUT_DIR="/wynton/group/databases/1000G_2504_high_coverage/working/20201028_3202_phased"
SAMPLES_FILE="/wynton/home/capra/gramey02/dnd_project/data/1KG_sample_names/unrelated_2504_samples.txt"
POP_FILE="/wynton/home/capra/gramey02/dnd_project/data/1KG_sample_names/unrelated_sample_superpop.txt"

BIALLELIC_DIR="/wynton/home/capra/gramey02/dnd_project/data/biallelic_vcfs_high_coverage"
AF_DIR="/wynton/home/capra/gramey02/dnd_project/data/biallelic_afs_high_coverage"

mkdir -p "$BIALLELIC_DIR"
mkdir -p "$AF_DIR"

# Map SGE_TASK_ID (1-23) -> chromosome label (1-22, X)
#CHROMS=({1..22} X)
CHR="X" #"${CHROMS[$((SGE_TASK_ID - 1))]}"

echo "=== Processing chr${CHR} ==="

INPUT_VCF="${INPUT_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.eagle2-phased.v2.vcf.gz"
# INPUT_VCF="${INPUT_DIR}/CCDG_14151_B01_GRM_WGS_2020-08-05_chr${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz"
BIALLELIC_VCF="${BIALLELIC_DIR}/TGP_chr${CHR}_biallelicSNPs.vcf.gz"
AF_FILE="${AF_DIR}/TGP_chr${CHR}_afs.txt"

if [[ ! -f "$INPUT_VCF" ]]; then
    echo "  ERROR: input file not found: $INPUT_VCF"
    exit 1
fi
if [[ ! -f "$SAMPLES_FILE" ]]; then
    echo "  ERROR: samples file not found: $SAMPLES_FILE"
    exit 1
fi
if [[ ! -f "$POP_FILE" ]]; then
    echo "  ERROR: population file not found: $POP_FILE"
    exit 1
fi

# Step 1: filter to biallelic SNPs
# Step 2: subset to unrelated samples in the same pass (-S)
# Step 3: recompute AC/AN/AF overall AND per-superpopulation (AC_EAS, AF_EAS, etc.
#         generated automatically from -S "$POP_FILE")
# Step 4: rename chr${CHR} -> ${CHR} in both header and data lines via --rename-chrs
bcftools view -m2 -M2 -v snps -S "$SAMPLES_FILE" --force-samples "$INPUT_VCF" \
    | bcftools +fill-tags -Ou -- -S "$POP_FILE" -t AC,AN,AF \
    | bcftools annotate --rename-chrs <(printf "chr%s\t%s\n" "$CHR" "$CHR") -Oz -o "$BIALLELIC_VCF"
bcftools index -t "$BIALLELIC_VCF"

# Step 6: extract fields (no header, since bcftools query never prints one)
bcftools query \
    -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AC\t%INFO/AN\t%INFO/AF\t%INFO/AF_EAS\t%INFO/AF_EUR\t%INFO/AF_AFR\t%INFO/AF_AMR\t%INFO/AF_SAS\n' \
    "$BIALLELIC_VCF" \
    > "$AF_FILE"

echo "  Wrote: $BIALLELIC_VCF"
echo "  Wrote: $AF_FILE"
echo "=== Done ==="