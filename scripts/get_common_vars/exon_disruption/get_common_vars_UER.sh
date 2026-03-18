#!/bin/bash
#$ -N get_common_vars_UER
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o /wynton/protected/home/capra/gramey02/ConklinCollab/scripts/out/get_common_vars_UER2.out
#$ -e /wynton/protected/home/capra/gramey02/ConklinCollab/scripts/err/get_common_vars_UER2.err

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
project_root="$(cd "$script_dir/../../.." && pwd)"

param_file="$2"
source "$param_file"
exon_file="$4"
output_dir="$1"
total_num_chroms="$3"

if [[ "$OUTPUT_DIR" = /* ]]; then
  resolved_output_base="$OUTPUT_DIR"
else
  resolved_output_base="$project_root/$OUTPUT_DIR"
fi

# col=$(head -1 $exon_file | tr ',' '\n' | grep -nx "chromosome_name" | cut -d: -f1)
# cut -d',' -f"$col" $exon_file | tail -n +2 | sort -u > $output_dir"/chromosomes/chrom_set.txt"
chrom_set="$resolved_output_base$RUN_NAME/chromosomes/chrom_set.txt"
#chrom_set=$CHROM_SET
chrom=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' "$chrom_set")
echo "$chrom"
# get the file that contains the biallelic snps on the current chromosome
af_file="/wynton/protected/home/capra/gramey02/ConklinCollab/data/dHS_and_related_GeneSets/Original_GeneSets/2025_04_22/chrom_separated_genes/vcfs/biallelic_afs/TGP_chr${chrom}_afs.txt"

# run the indel script
script="$script_dir/get_common_vars.py"
echo "Started identifying common coding vars..."
python3 "$script" --af_limit "$AF_LIMIT" --chrom "$chrom" --af_file "$af_file" --exon_file "$exon_file" --output_dir "$output_dir" --total_num_chroms "$total_num_chroms"
echo "Finished running chromosme ${chrom}."
