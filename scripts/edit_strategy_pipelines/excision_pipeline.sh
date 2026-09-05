#!/bin/bash
#$ -N excision_pipeline
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# parse input arguments
param_file="$2"
source "$param_file"
output_dir="$1"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

# set exon file
exon_file="$EXON_FILE_FOR_ANALYSIS"

module load CBI bcftools

# script to generate appropriate coordinates to look for excision variants in
find_variant_region="$project_root/scripts/format_variants/find_variant_region.py"
python3 "$find_variant_region" --exon_file "$project_root/$exon_file" --output_dir "$output_dir/CommonVars" --upstream_excision_dist "$UPSTREAM_EXCISION_DIST" --downstream_excision_dist "$DOWNSTREAM_EXCISION_DIST" --nearest_gene_file "$GENE_BODY_FILE" --filter_out_nearby_genes "$FILTER_OUT_NEARBY_GENES"

# script to generate list of variants in excision windows to run through excavate
find_commonVars="$project_root/scripts/get_common_vars/excision/find_excision_commonVars.sh"
gene_info="$output_dir/CommonVars/gene_filtered_excision_coords.txt"
echo "Started extracting common var info..."
qsub -sync y -l mem_free=1G -l h_rt=00:45:00 -o "$project_root/logs/out/find_excision_commonVars.out" -e "$project_root/logs/err/find_excision_commonVars.err" "$find_commonVars" "$output_dir/CommonVars" "$param_file" "$gene_info" "$exon_file"
echo "Finished extracting common variant info."

# get number of genes that have common vars in excision_window --> make sure there are 2 or more, since that's the minimum required for excision
common_var_genes="$output_dir/CommonVars/CommonVars_VarNumOver2_summary_noIDX.txt"
awk -F'\t' '$2+0 >= 2' "$output_dir/CommonVars/CommonVars_ALL_summary_noIDX.txt" > "$common_var_genes"
num_common_var_genes=$(wc -l < "$common_var_genes") # get the number of genes that have common vars in them (from summary file)

# script to determine viable pairs of snps, to narrow down the number of snps for future steps
filter_excision_snps="$project_root/scripts/format_variants/filter_excision_snps.sh"
cv_dict_filepath="$output_dir/CommonVars/CommonVars_ALL_dict.pkl"
echo "Filtering snps based on those that encompass exons across transcripts..."
qsub -t 1-"$num_common_var_genes" -sync y -l mem_free=10G -l h_rt=03:00:00 -o "$project_root/logs/out/filter_excision_snps.out" -e "$project_root/logs/err/filter_excision_snps.err" "$filter_excision_snps" "$output_dir/CommonVars/refined_common_vars" "$param_file" "$cv_dict_filepath" "$exon_file" "$common_var_genes" "$output_dir/CommonVars/valid_snp_pairs"
echo "Finished filtering excision snps."

# script to generate text files that we'll use to filter vcfs (which we'll then input into EXCAVATE)
generate_variant_textFiles="$project_root/scripts/format_variants/generate_variant_textFiles.py"
echo "Started generating common var loc files..."
#cv_dict_filepath="$output_dir/CommonVars/CommonVars_ALL_dict.pkl"
cv_dict_filepath="$output_dir/CommonVars/refined_common_vars"
python3 "$generate_variant_textFiles" --cv_dict_filepath "$cv_dict_filepath" --exon_file "$project_root/$exon_file" --output_dir "$output_dir" --af_file_dir "$project_root/$AF_FILE_DIR" --edit_strat "excision"
num_common_var_genes=$(ls -1 "$output_dir/excavate/CommonVar_locs" | wc -l)
echo "Finished generating common var loc files."

# script to filter vcfs accordingly
excavate_vcf_creation="$project_root/scripts/format_variants/generate_filtered_vcfs.sh" # first we need to generate vcf.gz files for the genes that have variants in their ubiquitous regions
input_metadata="$output_dir/excavate/input_metadata/excavate_run_metadata.txt"
echo "Started creating vcf files for excavate input..."
qsub -t 1-"$num_common_var_genes" -l mem_free=2G -l h_rt=01:00:00 -sync y -o "$project_root/logs/out/filt_vcfs_excision.out" -e "$project_root/logs/err/filt_vcfs_excision.err" "$excavate_vcf_creation" "$output_dir" "$param_file" "$input_metadata"
echo "Finished creating excavate inputs."

# split genes with too many variants into batches so EXCAVATE doesn't freeze/time out on them,
# run those batches through EXCAVATE separately, then merge each gene's batch outputs back together
batch_script="$project_root/scripts/format_variants/batch_excess_vars.py"
run_excavate_script="$project_root/scripts/excavate/run_excavate.sh"
merge_script="$project_root/scripts/format_variants/merge_batched_excavate_outputs.py"
unbatched_metadata="$output_dir/excavate/input_metadata/excavate_run_metadata_unbatched.txt"
batch_metadata="$output_dir/excavate/input_metadata/excavate_batch_metadata.txt"

if [[ "$DOWNSAMPLE" == "True" ]]; then
    echo "Started splitting high-variant genes into batches..."
    python3 "$batch_script" --output_dir "$output_dir" --ds_threshold "$DS_THRESHOLD"
    echo "Finished splitting high-variant genes into batches."

    if [[ -s "$batch_metadata" ]]; then
        num_batches=$(wc -l < "$batch_metadata")

        echo "Started creating vcf files for batched excavate input..."
        qsub -t 1-"$num_batches" -l mem_free=2G -l h_rt=01:00:00 -sync y -o "$project_root/logs/out/filt_vcfs_excision_batched.out" -e "$project_root/logs/err/filt_vcfs_excision_batched.err" "$excavate_vcf_creation" "$output_dir" "$param_file" "$batch_metadata"
        echo "Finished creating batched excavate inputs."

        echo "Started running EXCAVATE on batches..."
        qsub -t 1-"$num_batches" -l mem_free=4G -l h_rt=00:45:00 -sync y -o "$project_root/logs/out/excavate_excision_batched.out" -e "$project_root/logs/err/excavate_excision_batched.err" "$run_excavate_script" "$output_dir" "$param_file" "$batch_metadata"
        echo "Finished running EXCAVATE on batches."

        echo "Started merging batched EXCAVATE outputs back into per-gene folders..."
        python3 "$merge_script" --output_dir "$output_dir"
        echo "Finished merging batched EXCAVATE outputs."
    fi
    excavate_gene_info="$unbatched_metadata"
    num_excavate_genes=$(wc -l < "$unbatched_metadata")
else
    excavate_gene_info="$input_metadata"
    num_excavate_genes="$num_common_var_genes"
fi

# script to run excavate on the genes that didn't need batching
echo "Started running EXCAVATE..."
qsub -t 1-"$num_excavate_genes" -l mem_free=4G -l h_rt=00:45:00 -sync y -o "$project_root/logs/out/excavate_excision.out" -e "$project_root/logs/err/excavate_excision.err" "$run_excavate_script" "$output_dir" "$param_file" "$excavate_gene_info"
echo "Finished running EXCAVATE."

# generate text files for the valid guides so you can filter the individuals' vcfs accordingly
generate_guide_textFiles="$project_root/scripts/format_variants/generate_guide_textFiles.py"
guides_filepath="$output_dir/excavate/excavate_outputs"
echo "Started generating position files for valid guides..."
python3 "$generate_guide_textFiles" --guides_filepath "$guides_filepath" --exon_file "$project_root/$exon_file" --output_dir "$output_dir"
echo "Finished generating position files for valid guides."

# filter vcfs based on these viable EXCAVATE guides to see how many heterozygous individuals they will capture
guide_based_filtering="$project_root/scripts/format_variants/position_filtering.sh"
genes_w_guides="$output_dir/excavate/het_individuals/metadata/genes_w_valid_guides.txt"
num_genes_w_guides=$(awk -F'\t' '$1 != "" {n++} END{print n}' "$genes_w_guides")
echo "Started filtering vcfs based on valid guides..."
qsub -t 1-"$num_genes_w_guides" -l mem_free=3G -l h_rt=02:00:00 -sync y -o "$project_root/logs/out/guide_filtering_excision.out" -e "$project_root/logs/err/guide_filtering_excision.err" "$guide_based_filtering" "$output_dir" "$param_file" "$genes_w_guides"
echo "Finished filtering vcfs based on valid guides."

# script to calculate the number of heterozygous individuals for each gene, pre-excavate filtering
get_targeted_hets_prePAM="$project_root/scripts/get_hets/het_combos_prePAM.sh"
echo "Started calculating heterozygous individual numbers (prePAM filtering)..."
filtered_vcf_dir="$output_dir/excavate/input_vcfs"
gene_info="$output_dir/excavate/input_metadata/excavate_run_metadata.txt"
num_genes_prepam=$(awk -F'\t' '$1 != "" {n++} END{print n}' "$gene_info") # get number of prepam-filtered targetable genes from gene_info
valid_pairs_fp="$output_dir/CommonVars/valid_snp_pairs"
qsub -t 1-"$num_genes_prepam" -l mem_free=2G -l h_rt=01:00:00 -o "$project_root/logs/out/hets_excision_prepam_${RUN_NAME}.out" -e "$project_root/logs/err/hets_excision_prepam_${RUN_NAME}.err" "$get_targeted_hets_prePAM" "$output_dir/prePAM_hets" "$param_file" "$gene_info" "$filtered_vcf_dir" "$valid_pairs_fp"

# script to calculate number of heterozygous individuals that will be hit by pairs of snps
het_combos_script="$project_root/scripts/get_hets/het_combos.sh"
filtered_vcf_dir="$output_dir/excavate/Guide_filtered_vcfs"
valid_pairs_fp="$output_dir/CommonVars/valid_snp_pairs"
genes_w_guides="$output_dir/excavate/het_individuals/metadata/genes_w_valid_guides.txt"

# route genes with an outsized number of valid SNP pairs (e.g. DCC, RYR2) to a separate
# array job with a much larger mem_free/h_rt, since the greedy pair-selection algorithm
# in het_combos.py scales badly with pair count
split_script="$project_root/scripts/format_variants/split_genes_by_pair_count.py"
genes_w_guides_normal="$output_dir/excavate/het_individuals/metadata/genes_w_valid_guides_normal.txt"
genes_w_guides_large="$output_dir/excavate/het_individuals/metadata/genes_w_valid_guides_large.txt"
python3 "$split_script" --genes_w_guides "$genes_w_guides" --valid_pairs_fp "$valid_pairs_fp" \
    --threshold "$LARGE_GENE_PAIR_THRESHOLD" --normal_out "$genes_w_guides_normal" --large_out "$genes_w_guides_large"

echo "Started capturing heterozygote excision information..."

num_genes_normal=$(awk -F'\t' '$1 != "" {n++} END{print n}' "$genes_w_guides_normal")
if [[ "$num_genes_normal" -gt 0 ]]; then
    qsub -t 1-"$num_genes_normal" -l mem_free=2G -l h_rt=01:00:00 -sync y -o "$project_root/logs/out/hets_excision.out" -e "$project_root/logs/err/hets_excision.err" "$het_combos_script" "$output_dir/excavate/het_individuals" "$param_file" "$genes_w_guides_normal" "$filtered_vcf_dir" "$exon_file" "$valid_pairs_fp"
fi

num_genes_large=$(awk -F'\t' '$1 != "" {n++} END{print n}' "$genes_w_guides_large")
if [[ "$num_genes_large" -gt 0 ]]; then
    qsub -t 1-"$num_genes_large" -l mem_free=32G -l h_rt=08:00:00 -sync y -o "$project_root/logs/out/hets_excision_large.out" -e "$project_root/logs/err/hets_excision_large.err" "$het_combos_script" "$output_dir/excavate/het_individuals" "$param_file" "$genes_w_guides_large" "$filtered_vcf_dir" "$exon_file" "$valid_pairs_fp"
fi
echo "Finished capturing heterozygote excision information."
