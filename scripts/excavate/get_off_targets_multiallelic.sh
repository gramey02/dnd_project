#!/bin/bash
#$ -N get_off_targets_multiallelic
#$ -M Grace.Ramey@ucsf.edu
#$ -cwd

# Driver script that runs EXCAVATE with --off-targets across all editing strategies.
#
#   stage 1: non-excision strategies (acceptor_base_edits, CRISPRoff, donor_base_edits,
#            indels) run on their existing input_vcfs, since they carry few variants per
#            gene (median 1-2, max 19). All four are submitted at once WITHOUT -sync y so
#            they run concurrently with each other and with stages 2-3; the script
#            blocks on them at the very end instead.
#
#            NOTE: indel_pipeline.sh also calls run_excavate_w_offtargets.sh for indels.
#            That call has never produced output, but if it is ever re-enabled both would
#            write into indels/excavate/excavate_outputs_w_offtargets/.
#
#   stage 2: excision would be far too slow on its full variant set (135k variants across
#            537 genes; 171 genes over the batching threshold), so instead we build a
#            parallel results tree containing only the SNP loci that actually appear in
#            each gene's cross-strategy excision guide CSV (~26 loci/gene, max 47).
#
#   stage 3: run EXCAVATE with off-targets on those filtered excision VCFs.
#
# Nothing here modifies the existing excavate_outputs/ directories: every strategy's
# off-target results land in its own excavate/excavate_outputs_w_offtargets/ folder,
# and the excision run writes to excision/EXCISION_SELECTED_DIRNAME/ (below), a
# subfolder of excision/ that sits alongside the untouched excision/excavate/ tree.
#
# usage: qsub|bash get_off_targets_multiallelic.sh <run_dir> <param_file>
#   run_dir     e.g. /wynton/home/capra/gramey02/dnd_project/results/RUN_MULTIALLELIC
#   param_file  e.g. <run_dir>/PARAMS/params.txt
#
# set DRY_RUN=1 to print the qsub commands and build nothing, e.g.
#   DRY_RUN=1 bash get_off_targets_multiallelic.sh <run_dir> <param_file>

set -eo pipefail

# ---------------------------------------------------------------- input args
run_dir="$1"
param_file="$2"

if [[ -z "$run_dir" || -z "$param_file" ]]; then
    echo "usage: $0 <run_dir> <param_file>" >&2
    exit 1
fi
if [[ ! -f "$param_file" ]]; then
    echo "ERROR: param file not found: $param_file" >&2
    exit 1
fi

source "$param_file"
project_root="$PROJECT_ROOT"
script_dir="$project_root/scripts"

module load CBI miniforge3 bcftools

run_excavate_offtargets="$script_dir/excavate/run_excavate_w_offtargets.sh"
generate_filtered_vcfs="$script_dir/format_variants/generate_filtered_vcfs.sh"

# strategies that run on their full variant set
NON_EXCISION_STRATS=("acceptor_base_edits" "CRISPRoff" "donor_base_edits" "indels")

# where the SNP-restricted excision inputs/outputs get built: a subfolder of excision/,
# alongside the existing excision/excavate/ tree rather than inside it, so the original
# excision/excavate/ contents are never touched.
EXCISION_SELECTED_DIRNAME="excision_selected_snps_w_offtargets"
excision_dir="$run_dir/excision"
selected_dir="$excision_dir/$EXCISION_SELECTED_DIRNAME"

# same path relative to project_root, used only to fill in the trailing (unused) path
# columns of the generated metadata file. Handles run_dir being given either absolute
# or already relative to project_root.
selected_dir_rel="${selected_dir#"$project_root"/}"

# per-gene guide CSVs that define which excision SNP loci we care about
excision_guides_dir="$run_dir/summary_files/cross_strat_gRNAs/excision_guides/results"

# resource asks: off-target analysis scans the whole-genome fasta per guide set, so
# memory is driven by the genome load rather than by gene size. Matches the values
# indel_pipeline.sh already uses for run_excavate_w_offtargets.sh.
OFFTARGET_MEM="20G"
OFFTARGET_RT="08:00:00"
VCF_MEM="2G"
VCF_RT="01:00:00"

# ------------------------------------------------------------------- helpers

# job ids of the stage 1 arrays, which are submitted without -sync y so they run
# concurrently with each other and with stages 2-3. Waited on at the end.
stage1_jobids=()

submit() {
    # submit <mode> <description> <logname> <num_tasks> <mem> <runtime> -- <script> <args...>
    #   mode "sync"  - block until the array finishes (later stages depend on it)
    #   mode "async" - return as soon as it is queued, recording the job id
    local mode="$1" desc="$2" logname="$3" ntasks="$4" mem="$5" rt="$6"
    shift 6
    [[ "$1" == "--" ]] && shift

    if [[ "$ntasks" -le 0 ]]; then
        echo "  nothing to submit for $desc - skipping."
        return 0
    fi

    local cmd=(qsub -t "1-${ntasks}" -l "mem_free=${mem}" -l "h_rt=${rt}")
    [[ "$mode" == "sync" ]] && cmd+=(-sync y)
    cmd+=(-o "$project_root/logs/out/${logname}.out"
          -e "$project_root/logs/err/${logname}.err"
          "$@")

    if [[ -n "$DRY_RUN" ]]; then
        echo "  [DRY_RUN] ${cmd[*]}"
        return 0
    fi

    echo "  submitting $ntasks tasks for $desc ..."
    if [[ "$mode" == "sync" ]]; then
        "${cmd[@]}"
        echo "  finished $desc."
    else
        local qsub_out jobid
        qsub_out=$("${cmd[@]}")
        echo "    $qsub_out"
        # 'Your job-array 12345.1-136:1 ("name") has been submitted' -> 12345
        jobid=$(echo "$qsub_out" | awk '{print $3}' | cut -d. -f1)
        stage1_jobids+=("$jobid")
        echo "  queued $desc as job $jobid - not waiting."
    fi
}

# =============================================================== STAGE 1 ====
# non-excision strategies, run on their existing (unfiltered) input VCFs
echo "=== STAGE 1: off-targets for non-excision strategies ==="

for strat in "${NON_EXCISION_STRATS[@]}"; do
    strat_dir="$run_dir/$strat"
    strat_metadata="$strat_dir/excavate/input_metadata/excavate_run_metadata.txt"

    if [[ ! -s "$strat_metadata" ]]; then
        echo "  WARNING: no metadata at $strat_metadata - skipping $strat." >&2
        continue
    fi

    num_genes=$(wc -l < "$strat_metadata")
    echo "[$strat] $num_genes genes"
    submit async "$strat" "excavate_${strat}_w_offtargets" "$num_genes" \
        "$OFFTARGET_MEM" "$OFFTARGET_RT" -- \
        "$run_excavate_offtargets" "$strat_dir" "$param_file" "$strat_metadata"
done

echo "  all non-excision strategies queued; continuing to excision while they run."

# =============================================================== STAGE 2 ====
# build a SNP-restricted copy of the excision inputs
echo
echo "=== STAGE 2: building SNP-restricted excision inputs ==="

excision_metadata="$excision_dir/excavate/input_metadata/excavate_run_metadata.txt"
selected_metadata="$selected_dir/excavate/input_metadata/excavate_run_metadata.txt"

for sub in CommonVar_locs input_vcfs input_metadata excavate_outputs_w_offtargets; do
    mkdir -p "$selected_dir/excavate/$sub"
done

# Pull the SNP positions out of each gene's guide CSV and write them in the same
# two-column "<chrom>\t<pos>" format that generate_filtered_vcfs.sh expects, then
# emit a matching metadata row carrying the gene's FULL locus (unchanged) so guide
# coordinates stay comparable to the original excision run.
export GUIDES_DIR="$excision_guides_dir"
export SRC_METADATA="$excision_metadata"
export OUT_LOCS_DIR="$selected_dir/excavate/CommonVar_locs"
export OUT_METADATA="$selected_metadata"
export OUT_REL_PREFIX="$selected_dir_rel"

python3 <<'PYEOF'
import csv
import glob
import os

guides_dir = os.environ["GUIDES_DIR"]
src_metadata = os.environ["SRC_METADATA"]
out_locs_dir = os.environ["OUT_LOCS_DIR"]
out_metadata = os.environ["OUT_METADATA"]
rel_prefix = os.environ["OUT_REL_PREFIX"]

# gene -> full metadata row from the original (unbatched-locus) excision metadata
meta = {}
with open(src_metadata) as fh:
    for line in fh:
        if not line.strip():
            continue
        fields = line.rstrip("\n").split("\t")
        meta[fields[0]] = fields

csv_files = sorted(glob.glob(os.path.join(guides_dir, "*_excision_gRNAs.csv")))
if not csv_files:
    raise SystemExit(f"ERROR: no guide CSVs found in {guides_dir}")

written, skipped_no_meta, skipped_no_snps, total_snps = [], [], [], 0

for path in csv_files:
    gene = os.path.basename(path).replace("_excision_gRNAs.csv", "")

    if gene not in meta:
        skipped_no_meta.append(gene)
        continue

    # snp columns look like "<position>_<allele_index>", e.g. 70269677_0
    positions = set()
    with open(path) as fh:
        for row in csv.DictReader(fh):
            for col in ("snp1_allele", "snp2_allele"):
                value = (row.get(col) or "").strip()
                if value:
                    positions.add(int(value.split("_")[0]))

    if not positions:
        skipped_no_snps.append(gene)
        continue

    chrom = meta[gene][1]
    with open(os.path.join(out_locs_dir, f"{gene}_CommonVar_locs.txt"), "w") as out:
        for pos in sorted(positions):
            out.write(f"{chrom}\t{pos}\n")

    written.append(gene)
    total_snps += len(positions)

# rewrite the trailing path columns so the metadata points into the new tree
with open(out_metadata, "w") as out:
    for gene in written:
        fields = list(meta[gene])
        while len(fields) < 7:
            fields.append("")
        fields[5] = (
            f"{rel_prefix}/excavate/CommonVar_locs/{gene}_CommonVar_locs.txt"
        )
        fields[6] = (
            f"{rel_prefix}/excavate/input_vcfs/{gene}_CommonVar_filtered.vcf.gz"
        )
        out.write("\t".join(fields) + "\n")

print(f"  guide CSVs found:       {len(csv_files)}")
print(f"  genes written:          {len(written)}")
print(f"  total SNP loci kept:    {total_snps}")
if written:
    print(f"  mean SNP loci per gene: {total_snps / len(written):.1f}")
if skipped_no_meta:
    print(f"  SKIPPED (no excision metadata row): {len(skipped_no_meta)} -> "
          f"{', '.join(skipped_no_meta[:10])}")
if skipped_no_snps:
    print(f"  SKIPPED (no SNPs in guide CSV):     {len(skipped_no_snps)} -> "
          f"{', '.join(skipped_no_snps[:10])}")
PYEOF

if [[ ! -s "$selected_metadata" ]]; then
    echo "ERROR: no genes written to $selected_metadata - nothing to run for excision." >&2
    exit 1
fi

num_selected_genes=$(wc -l < "$selected_metadata")

# filter each gene's chromosome VCF down to just those SNP loci. Reuses
# generate_filtered_vcfs.sh unmodified: it reads CommonVar_locs/ and writes
# input_vcfs/ relative to whatever output_dir it is handed.
echo "[excision] creating filtered VCFs for $num_selected_genes genes"
submit sync "excision filtered VCFs" "filt_vcfs_excision_selected_snps" "$num_selected_genes" \
    "$VCF_MEM" "$VCF_RT" -- \
    "$generate_filtered_vcfs" "$selected_dir" "$param_file" "$selected_metadata"

# =============================================================== STAGE 3 ====
# off-targets on the SNP-restricted excision inputs
echo
echo "=== STAGE 3: off-targets for excision (SNP-restricted) ==="
echo "[excision] $num_selected_genes genes"
submit sync "excision off-targets" "excavate_excision_selected_snps_w_offtargets" \
    "$num_selected_genes" "$OFFTARGET_MEM" "$OFFTARGET_RT" -- \
    "$run_excavate_offtargets" "$selected_dir" "$param_file" "$selected_metadata"

# ------------------------------------------------------------------ barrier
# stages 2-3 are independent of stage 1, so the non-excision arrays have been
# running this whole time. Block on them before reporting DONE by submitting a
# trivial job that is held until every stage 1 job id has finished.
if [[ ${#stage1_jobids[@]} -gt 0 ]]; then
    echo
    echo "=== waiting on stage 1 jobs: ${stage1_jobids[*]} ==="
    hold_list=$(IFS=,; echo "${stage1_jobids[*]}")
    qsub -sync y -hold_jid "$hold_list" -b y -N off_targets_barrier \
        -l mem_free=1G -l h_rt=00:05:00 \
        -o /dev/null -e /dev/null /bin/true
    echo "  stage 1 complete."
fi

echo
echo "=== DONE ==="
for strat in "${NON_EXCISION_STRATS[@]}"; do
    echo "  $strat -> $run_dir/$strat/excavate/excavate_outputs_w_offtargets/"
done
echo "  excision      -> $selected_dir/excavate/excavate_outputs_w_offtargets/"
