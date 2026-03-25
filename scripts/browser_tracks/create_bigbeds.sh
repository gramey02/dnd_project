#!/bin/bash
#$ -N create_bigbeds
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ../../logs/out/create_bigbeds.out
#$ -e ../../logs/err/create_bigbeds.err

# Set directories
bt_dir="../../data/browser_tracks/"
results_dir="$bt_dir/bed_files"
# Chrom sizes file (REQUIRED for bigBed)
chrom_sizes="hg38.chrom.sizes"
# file with all targetable gene names
gene_file="../../results/summary_files/targetable_genes.txt"
# base directory containing per-gene files
pg_dir="$results_dir/per_gene_files"
# # isolate one gene at a time for an array job
# gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $gene_file)

# Public URLs reused across the UCSC hub text files
github_base_url="https://raw.githubusercontent.com/gramey02/DnD_TrackHubs_Public/refs/heads/main/bed_files"
cN8_url="https://www.dropbox.com/scl/fi/v0zhxb1nv0s3lbnvtsfdi/cN8-hNIL.vcf.gz?rlkey=5jvvckcl5irvgquergenz622v&st=f8k98eqe&dl=0"
cN8_tbi="https://www.dropbox.com/scl/fi/3vvzenwaeptn4bk5uq3p8/cN8-hNIL.vcf.gz.tbi?rlkey=cdjowthr6v2ipn6lb11ffvdu2&st=uimb4dce&dl=0"
KOLF2_url="https://www.dropbox.com/scl/fi/khggb81kz3boboa08i5ao/KOLF2-ARID2-A02.vcf.gz?rlkey=plw72thalgc18clhllj2qt9yu&st=razkh9zp&dl=0"
KOLF2_tbi="https://www.dropbox.com/scl/fi/wjnmh6ooqzsn41cej4vkf/KOLF2-ARID2-A02.vcf.gz.tbi?rlkey=8g5edij1dsgpdgb0pl6cos93k&st=9iglrkqt&dl=0"
WTB_url="https://www.dropbox.com/scl/fi/2convcch36hvpc6yrkck7/WTB_variants_PASS.vcf.gz?rlkey=xpzy5409qwppftdssnt5puytp&st=z9nxdvt2&dl=0"
WTB_tbi="https://www.dropbox.com/scl/fi/msnbr6rbgx6nii4tvqpcs/WTB_variants_PASS.vcf.gz.tbi?rlkey=qa64d30avsprcnkn3ija51kck&st=6x16hnr9&dl=0"
WTC_url="https://www.dropbox.com/scl/fi/xtbksv9x1ufdooditcyel/WTC_variants_PASS.vcf.gz?rlkey=8hvqxbo4mycyh146n2betc54k&st=johlzb9u&dl=0"
WTC_tbi="https://www.dropbox.com/scl/fi/dvn57oncyxyuelyw79nu2/WTC_variants_PASS.vcf.gz.tbi?rlkey=8ntna5anht7uju2505yoqnhqx&st=eh54x6ar&dl=0"

write_gene_hub_file() {
    local gene="$1"
    local out_path="$2"
    # Point this gene's hub file at its per-gene bigBed on GitHub.
    local bigbed_url="${github_base_url}/per_gene_files/${gene}/${gene}_ng.bb"

    # Write a single-track UCSC hub config for this gene.
    cat > "$out_path" <<EOF
hub Dominant & Dispensible Gene Editing Opportunities - ${gene}
shortLabel D&D Gene Editing - ${gene}
longLabel Common genetic variant hub for dominant and dispensible (D&D) gene editing opportunities (${gene})
useOneFile on
email Grace.Ramey@ucsf.edu

genome hg38

track ${gene} Common Variant Editing Targets
shortLabel ${gene} targets
longLabel Mutation-agnostic and allele-specific editing sites for ${gene}
visibility pack
type bigBed
bigDataUrl ${bigbed_url}
EOF
}

write_gene_hub_file_with_cell_lines() {
    local gene="$1"
    local out_path="$2"
    # Reuse the same per-gene bigBed URL, then add external VCF tracks.
    local bigbed_url="${github_base_url}/per_gene_files/${gene}/${gene}_ng.bb"

    # Write a gene hub config that includes the bigBed plus cell-line VCFs.
    cat > "$out_path" <<EOF
hub Dominant & Dispensible Gene Editing Opportunities - ${gene}
shortLabel D&D Gene Editing - ${gene}
longLabel Common genetic variant hub for dominant and dispensible (D&D) gene editing opportunities (${gene})
useOneFile on
email Grace.Ramey@ucsf.edu

genome hg38

track ${gene} Common Variant Editing Targets
shortLabel ${gene} targets
longLabel Mutation-agnostic and allele-specific editing sites for ${gene}
visibility pack
type bigBed
bigDataUrl ${bigbed_url}

track WTB_vcf
shortLabel WTB variants
longLabel Cell Line Variants: WTB
visibility pack
type vcfTabix
bigDataUrl ${WTB_url}

track WTC_vcf
shortLabel WTC variants
longLabel Cell Line Variants: WTC
visibility pack
type vcfTabix
bigDataUrl ${WTC_url}

track WTD_vcf
shortLabel WTD variants
longLabel Cell Line Variants: WTD
visibility pack
type vcfTabix
bigDataUrl ${cN8_url}

track KOLF2_vcf
shortLabel KOLF2 variants
longLabel Cell Line Variants: KOLF2
visibility pack
type vcfTabix
bigDataUrl ${KOLF2_url}
EOF
}

write_metadata_hub_file() {
    local out_path="$1"
    # This hub file points to the combined all-gene bigBed.
    local bigbed_url="${github_base_url}/metadata/DnD_gene_ng.bb"

    # Write the single-track UCSC hub config for the combined dataset.
    cat > "$out_path" <<EOF
hub Dominant & Dispensible Gene Editing Opportunities
shortLabel D&D Gene Editing Opportunities
longLabel Common genetic variant hub for dominant and dispensible (D&D) gene editing opportunities
useOneFile on
email Grace.Ramey@ucsf.edu

genome hg38

track CommonVar_EditingTargets
shortLabel D&D gene editing targets
longLabel Mutation-agnostic and allele-specific editing sites for D&D genes
visibility pack
type bigBed
bigDataUrl ${bigbed_url}
EOF
}

write_metadata_hub_file_with_cell_lines() {
    local out_path="$1"
    # The metadata hub uses the combined bigBed and indexed VCF tracks.
    local bigbed_url="${github_base_url}/metadata/DnD_gene_ng.bb"

    # Write the combined hub config plus cell-line VCF/index URLs.
    cat > "$out_path" <<EOF
hub Dominant & Dispensible Gene Editing Opportunities
shortLabel D&D Gene Editing Opportunities
longLabel Common genetic variant hub for dominant and dispensible (D&D) gene editing opportunities
useOneFile on
email Grace.Ramey@ucsf.edu

genome hg38

track CommonVar_EditingTargets
shortLabel D&D gene editing targets
longLabel Mutation-agnostic and allele-specific editing sites for D&D genes
visibility pack
type bigBed
bigDataUrl ${bigbed_url}

track WTB_vcf
shortLabel WTB variants
longLabel Cell Line Variants: WTB
visibility pack
type vcfTabix
bigDataUrl ${WTB_url}
bigDataIndex ${WTB_tbi}

track WTC_vcf
shortLabel WTC variants
longLabel Cell Line Variants: WTC
visibility pack
type vcfTabix
bigDataUrl ${WTC_url}
bigDataIndex ${WTC_tbi}

track WTD_vcf
shortLabel WTD variants
longLabel Cell Line Variants: WTD
visibility pack
type vcfTabix
bigDataUrl ${cN8_url}
bigDataIndex ${cN8_tbi}

track KOLF2_vcf
shortLabel KOLF2 variants
longLabel Cell Line Variants: KOLF2
visibility pack
type vcfTabix
bigDataUrl ${KOLF2_url}
bigDataIndex ${KOLF2_tbi}
EOF
}

# activate conda environment where bigBed conversion package is loaded
module load CBI miniforge3
conda activate dnscripts

# Create temp file for combined bed
combined_bed="$results_dir/metadata/DnD_gene.combined.bed"
> "$combined_bed"  # empty file

# Loop through genes to create per-gene bigBeds, hub text files, and the combined BED
while read gene; do
    gene_dir="$pg_dir/$gene"
    bed_file="$gene_dir/${gene}_snp_track_ng.bed"
    sorted_bed="$gene_dir/${gene}_snp_track_ng.sorted.bed"
    bigbed_file="$gene_dir/${gene}_ng.bb"
    as_file="$gene_dir/${gene}_ng.as"

    if [[ -f "$bed_file" ]]; then
        echo "Processing $gene"

        # Sort BED for bigBed
        sort -k1,1 -k2,2n "$bed_file" > "$sorted_bed"

        # Convert to bigBed
        bedToBigBed -type=bed9+15 -tab -as=$as_file $sorted_bed $chrom_sizes $bigbed_file

        # Append to combined BED
        cat "$sorted_bed" >> "$combined_bed"

        # Create the plain hub file and the cell-line-augmented hub file for this gene.
        write_gene_hub_file "$gene" "$gene_dir/${gene}_hub_file_ng.txt"
        write_gene_hub_file_with_cell_lines "$gene" "$gene_dir/${gene}_hub_file_wCellLines_ng.txt"

        rm "$sorted_bed"
    else
        echo "Missing BED file for $gene"
    fi
done < "$gene_file"

# Make the final combined bigBed
first_gene=$(head -n1 "$gene_file")
gene_dir="$pg_dir/$first_gene"
cp "$gene_dir/${first_gene}_ng.as" "$results_dir/metadata/"
as_file="$results_dir/metadata/${first_gene}_ng.as"
sorted_combined="$results_dir/metadata/DnD_gene_ng.sorted.bed"
final_bigbed="$results_dir/metadata/DnD_gene_ng.bb"

sort -k1,1 -k2,2n "$combined_bed" > "$sorted_combined"
bedToBigBed -type=bed9+15 -tab -as=$as_file $sorted_combined $chrom_sizes $final_bigbed

# Create the two top-level hub files for the combined all-gene tracks.
write_metadata_hub_file "$results_dir/metadata/All_DnD_genes_hub_file_ng.txt"
write_metadata_hub_file_with_cell_lines "$results_dir/metadata/All_DnD_genes_hub_file_wCellLines_ng.txt"

rm "$sorted_combined" "$combined_bed"

echo "Finished converting beds to bigBeds."
