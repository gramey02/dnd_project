#!/bin/bash
#$ -N create_bigbeds
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd
#$ -o ../../logs/out/create_bigbeds.out
#$ -e ../../logs/err/create_bigbeds.err

# Set directories
# parse input args
param_file="$1"
gene_file="$2"
output_dir="$3"
# load variables from input args
source $param_file
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"
run_name=$RUN_NAME

# set overall browser track (bt) directory
bt_dir="$output_dir/summary_files/browser_tracks"
# Chrom sizes file (REQUIRED for bigBed)
chrom_sizes="$project_root/data/chrom_sizes/hg38.chrom.sizes"
# base directory containing per-gene files
pg_dir="$bt_dir/per_gene_files"
# # isolate one gene at a time for an array job
# gene=$(awk -v row=$SGE_TASK_ID 'NR == row {print $1}' $gene_file)

# Public URLs reused across the UCSC hub text files
github_base_url="https://raw.githubusercontent.com/gramey02/DnD_TrackHubs_Public/refs/heads/main/bed_files"
cN8_url="https://raw.githubusercontent.com/gramey02/DnD_TrackHubs_Public/refs/heads/main/all_genes_w_filtering/WTD_het.bb" # "https://www.dropbox.com/scl/fi/v0zhxb1nv0s3lbnvtsfdi/cN8-hNIL.vcf.gz?rlkey=5jvvckcl5irvgquergenz622v&st=f8k98eqe&dl=0"
# cN8_tbi="https://raw.githubusercontent.com/gramey02/DnD_TrackHubs_Public/refs/heads/main/all_genes_w_filtering/WTD_het.bb" # "https://www.dropbox.com/scl/fi/3vvzenwaeptn4bk5uq3p8/cN8-hNIL.vcf.gz.tbi?rlkey=cdjowthr6v2ipn6lb11ffvdu2&st=uimb4dce&dl=0"
KOLF2_url="https://raw.githubusercontent.com/gramey02/DnD_TrackHubs_Public/refs/heads/main/all_genes_w_filtering/KOLF2_het.bb" # "https://www.dropbox.com/scl/fi/khggb81kz3boboa08i5ao/KOLF2-ARID2-A02.vcf.gz?rlkey=plw72thalgc18clhllj2qt9yu&st=razkh9zp&dl=0"
# KOLF2_tbi="https://www.dropbox.com/scl/fi/wjnmh6ooqzsn41cej4vkf/KOLF2-ARID2-A02.vcf.gz.tbi?rlkey=8g5edij1dsgpdgb0pl6cos93k&st=9iglrkqt&dl=0"
WTB_url="https://raw.githubusercontent.com/gramey02/DnD_TrackHubs_Public/refs/heads/main/all_genes_w_filtering/WTB_het.bb" # "https://www.dropbox.com/scl/fi/2convcch36hvpc6yrkck7/WTB_variants_PASS.vcf.gz?rlkey=xpzy5409qwppftdssnt5puytp&st=z9nxdvt2&dl=0"
# WTB_tbi="https://www.dropbox.com/scl/fi/msnbr6rbgx6nii4tvqpcs/WTB_variants_PASS.vcf.gz.tbi?rlkey=qa64d30avsprcnkn3ija51kck&st=6x16hnr9&dl=0"
WTC_url="https://raw.githubusercontent.com/gramey02/DnD_TrackHubs_Public/refs/heads/main/all_genes_w_filtering/WTC_het.bb" # "https://www.dropbox.com/scl/fi/xtbksv9x1ufdooditcyel/WTC_variants_PASS.vcf.gz?rlkey=8hvqxbo4mycyh146n2betc54k&st=johlzb9u&dl=0"
# WTC_tbi="https://www.dropbox.com/scl/fi/dvn57oncyxyuelyw79nu2/WTC_variants_PASS.vcf.gz.tbi?rlkey=8ntna5anht7uju2505yoqnhqx&st=eh54x6ar&dl=0"

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

track ${gene}_Common_Variant_Editing_Targets
shortLabel ${gene} targets
longLabel Mutation-agnostic and allele-specific editing sites for ${gene}
visibility pack
type bigBed
filterValues.editStrats exon disruption,epigenetic silencing,splice site disruption,excision
filterType.editStrats multipleListOr
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

track ${gene}_Common_Variant_Editing_Targets
shortLabel ${gene} targets
longLabel Mutation-agnostic and allele-specific editing sites for ${gene}
visibility pack
type bigBed
filterValues.editStrats exon disruption,epigenetic silencing,splice site disruption,excision
filterType.editStrats multipleListOr
bigDataUrl ${bigbed_url}

track WTB_bb
shortLabel WTB het vars
longLabel WTB Cell Line Heterozygous Variants
visibility pack
type bigBed
filterValues.snv_qual HIGH,LOW
filterType.snv_qual multipleListOr
bigDataUrl ${WTB_url}

track WTC_bb
shortLabel WTC het vars
longLabel WTC Cell Line Heterozygous Variants
visibility pack
type bigBed
filterValues.snv_qual HIGH,LOW
filterType.snv_qual multipleListOr
bigDataUrl ${WTC_url}

track WTD_bb
shortLabel WTD het vars
longLabel WTD Cell Line Heterozygous Variants
visibility pack
type bigBed
filterValues.snv_qual HIGH,LOW
filterType.snv_qual multipleListOr
bigDataUrl ${cN8_url}

track KOLF2_bb
shortLabel KOLF2 het vars
longLabel KOLF2 Cell Line Heterozygous Variants
visibility pack
type bigBed
filterValues.snv_qual HIGH,LOW
filterType.snv_qual multipleListOr
bigDataUrl ${KOLF2_url}
EOF
}

write_metadata_hub_file() {
    local out_path="$1"
    # This hub file points to the combined all-gene bigBed.
    local bigbed_url="${github_base_url}/all_genes_w_filtering/DnD_gene_ng.bb"

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
filterValues.editStrats exon disruption,epigenetic silencing,splice site disruption,excision
filterType.editStrats multipleListOr
bigDataUrl ${bigbed_url}
EOF
}

write_metadata_hub_file_with_cell_lines() {
    local out_path="$1"
    # The metadata hub uses the combined bigBed and indexed VCF tracks.
    local bigbed_url="${github_base_url}/all_genes_w_filtering/DnD_gene_ng.bb"

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
filterValues.editStrats exon disruption,epigenetic silencing,splice site disruption,excision
filterType.editStrats multipleListOr
bigDataUrl ${bigbed_url}

track WTB_bb
shortLabel WTB het vars
longLabel WTB Cell Line Heterozygous Variants
visibility pack
type bigBed
filterValues.snv_qual HIGH,LOW
filterType.snv_qual multipleListOr
bigDataUrl ${WTB_url}

track WTC_bb
shortLabel WTC het vars
longLabel WTC Cell Line Heterozygous Variants
visibility pack
type bigBed
filterValues.snv_qual HIGH,LOW
filterType.snv_qual multipleListOr
bigDataUrl ${WTC_url}

track WTD_bb
shortLabel WTD het vars
longLabel WTD Cell Line Heterozygous Variants
visibility pack
type bigBed
filterValues.snv_qual HIGH,LOW
filterType.snv_qual multipleListOr
bigDataUrl ${cN8_url}

track KOLF2_bb
shortLabel KOLF2 het vars
longLabel KOLF2 Cell Line Heterozygous Variants
visibility pack
type bigBed
filterValues.snv_qual HIGH,LOW
filterType.snv_qual multipleListOr
bigDataUrl ${KOLF2_url}
EOF
}

# activate conda environment where bigBed conversion package is loaded
module load CBI miniforge3
conda activate merged_env # dnscripts

# Locate combined bed file
combined_bed="$bt_dir/all_genes_together/all_genes.bed"

# Loop through genes to create per-gene bigBeds, hub text files, and the combined BED
while read gene; do
    gene_dir="$pg_dir/$gene"
    bed_file="$gene_dir/${gene}_snp_track_ng.bed"
    sorted_bed="$gene_dir/${gene}_snp_track_ng.sorted.bed"
    bigbed_file="$gene_dir/${gene}_ng.bb"
    as_file="$bt_dir/metadata/bed_col_descriptors.as"

    if [[ -f "$bed_file" ]]; then
        echo "Processing $gene"

        # Sort BED for bigBed
        sort -k1,1 -k2,2n "$bed_file" > "$sorted_bed"

        # Convert to bigBed
        bedToBigBed -type=bed9+15 -tab -as=$as_file $sorted_bed $chrom_sizes $bigbed_file

        # # Append to combined BED
        # cat "$sorted_bed" >> "$combined_bed"

        # Create the plain hub file and the cell-line-augmented hub file for this gene.
        write_gene_hub_file "$gene" "$gene_dir/${gene}_hub_file_ng.txt"
        write_gene_hub_file_with_cell_lines "$gene" "$gene_dir/${gene}_hub_file_wCellLines_ng.txt"

        rm "$sorted_bed"
    else
        echo "Missing BED file for $gene"
    fi
done < "$gene_file"

# Make the final combined bigBed
as_file="$bt_dir/metadata/bed_col_descriptors.as"
sorted_combined="$bt_dir/all_genes_together/all_genes_ng.sorted.bed"
final_bigbed="$bt_dir/all_genes_together/DnD_gene_ng.bb"

sort -k1,1 -k2,2n "$combined_bed" > "$sorted_combined"
bedToBigBed -type=bed9+15 -tab -as=$as_file $sorted_combined $chrom_sizes $final_bigbed

# Create the two top-level hub files for the combined all-gene tracks.
write_metadata_hub_file "$bt_dir/all_genes_together/All_DnD_genes_hub_file_ng.txt"
write_metadata_hub_file_with_cell_lines "$bt_dir/all_genes_together/All_DnD_genes_hub_file_wCellLines_ng.txt"

echo "Finished converting beds to bigBeds."
