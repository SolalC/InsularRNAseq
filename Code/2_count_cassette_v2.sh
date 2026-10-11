#!/bin/bash
#SBATCH --job-name=hfd_count_v2
#SBATCH --partition=general
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-12:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/countv2_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/countv2_%A_%a.err

# Step 2 (v2 reference): cellranger count, one array task per mouse.
#
# Identical to 2_count.sh except for the reference and the output root. Runs
# against GRCm39_cassette_v2, whose transgene contigs are cassette-derived: the
# v1 'shared' 3' contig was cut from p24 and matches only 395/799 bp of p131 and
# 589/799 of p140, so the 3' end of the tdTomato and EBFP2 constructs was absent
# from the reference. See Results/2_4_cassette_report/FINDINGS.md.
#
# Output goes to Data/Counts_v2 so Data/Counts (the v1 run) stays intact for
# comparison. The transgene feature names change: EGFP/tdTomato/EBFP2 are now
# whole-cassette contigs and the shared element is 'AAV_shared', not
# 'AAV_WPRE_SV40pA_shared' - downstream scripts keyed to the old name need updating.
#
# The three sequencing pools are the SAME libraries re-sequenced, so they are
# lanes of one library, not replicates: every pool contributes one lane
# (L001/L002/L003) under a single --sample prefix, and all three fastq
# directories are passed together to one cellranger run.
#
# Sample list is driven by Data/sample_metadata.csv (column fastq_prefix,
# rows with fastq_ok==True) so there is one source of truth for sample names.

set -euo pipefail

PROJ=/scratch/user/uqschauq/Insular
CELLRANGER=$PROJ/Soft/cellranger-10.1.0/bin/cellranger
REF=$PROJ/Data/Reference/GRCm39_cassette_v2
FASTQROOT=$PROJ/Data
MANIFEST=$PROJ/Data/sample_metadata.csv
OUTROOT=$PROJ/Data/Counts_v2

# snRNA-seq: intronic reads are most of the signal, so they are counted (this is
# the 10.1 default, set explicitly for the record).
INCLUDE_INTRONS=true
# BAM kept: the transgene calls are presence/absence at a 1-UMI threshold, and
# checking a suspect call (ambient? index hop? the 98 bp EGFP/EBFP2 identical
# run?) means going back to the reads. Set false to save ~40-60 GB per sample.
CREATE_BAM=true

mkdir -p "$OUTROOT" "$PROJ/Code/logs"

echo "[$(date)] host=$(hostname) job=${SLURM_ARRAY_JOB_ID:-none} task=${SLURM_ARRAY_TASK_ID:-none}"

# ---- inputs -----------------------------------------------------------------
[[ -x $CELLRANGER ]] || { echo "FATAL: cellranger not executable at $CELLRANGER"; exit 1; }
[[ -s $REF/reference.json ]] || { echo "FATAL: no built reference at $REF"; exit 1; }
[[ -s $MANIFEST ]] || { echo "FATAL: missing manifest $MANIFEST"; exit 1; }

# ---- which sample is this task ---------------------------------------------
IDX=${SLURM_ARRAY_TASK_ID:-1}
# the manifest has CRLF line endings and the columns are looked up by name, so
# neither a stray \r nor a reordered column can silently change the sample list
mapfile -t SAMPLES < <(awk -F, '
    { sub(/\r$/, "") }
    NR==1 { for (i=1; i<=NF; i++) col[$i]=i;
            if (!("fastq_prefix" in col) || !("fastq_ok" in col)) {
                print "FATAL: manifest lacks fastq_prefix/fastq_ok columns" > "/dev/stderr"; exit 1 }
            next }
    $(col["fastq_ok"]) == "True" { print $(col["fastq_prefix"]) }' "$MANIFEST")
[[ ${#SAMPLES[@]} -eq 12 ]] || { echo "FATAL: expected 12 usable samples in manifest, got ${#SAMPLES[@]}"; exit 1; }
[[ $IDX -ge 1 && $IDX -le ${#SAMPLES[@]} ]] || { echo "FATAL: array index $IDX outside 1..${#SAMPLES[@]}"; exit 1; }
SAMPLE=${SAMPLES[$((IDX-1))]}
OUT=$OUTROOT/$SAMPLE
echo "sample     : $SAMPLE  (manifest row $IDX)"

# ---- locate this sample's fastq directories --------------------------------
# Directory naming differs between pools (Ech1_FSD-ds.*, Ech1FSD_1-ds.*,
# Ech1FSD_2-ds.*), so find them by the fastq filenames they contain instead.
mapfile -t FQDIRS < <(find "$FASTQROOT" -maxdepth 2 -name "${SAMPLE}_S*_L00?_R1_001.fastq.gz" \
                        -printf '%h\n' | sort -u)
if [[ ${#FQDIRS[@]} -ne 3 ]]; then
    echo "FATAL: expected 3 fastq directories for $SAMPLE, found ${#FQDIRS[@]}:"
    printf '  %s\n' "${FQDIRS[@]}"
    exit 1
fi

# every lane must be distinct, and R1/R2 must both be present in each directory
LANES=()
for d in "${FQDIRS[@]}"; do
    for r in R1 R2; do
        n=$(find "$d" -maxdepth 1 -name "${SAMPLE}_S*_L00?_${r}_001.fastq.gz" | wc -l)
        [[ $n -eq 1 ]] || { echo "FATAL: $d has $n ${r} files for $SAMPLE, expected 1"; exit 1; }
    done
    LANES+=( "$(find "$d" -maxdepth 1 -name "${SAMPLE}_S*_R1_001.fastq.gz" -printf '%f\n' | grep -oE 'L00[0-9]')" )
done
UNIQ_LANES=$(printf '%s\n' "${LANES[@]}" | sort -u | wc -l)
[[ $UNIQ_LANES -eq 3 ]] || { echo "FATAL: lanes not distinct across pools for $SAMPLE: ${LANES[*]}"; exit 1; }
echo "lanes      : $(printf '%s ' "${LANES[@]}")"
printf 'fastq dir  : %s\n' "${FQDIRS[@]}"

FASTQS=$(IFS=,; echo "${FQDIRS[*]}")

# ---- refuse to clobber a previous run --------------------------------------
if [[ -e $OUT ]]; then
    if [[ -s $OUT/outs/metrics_summary.csv ]]; then
        echo "SKIP: $SAMPLE already has completed outs at $OUT — nothing to do."
        exit 0
    fi
    if [[ ${FORCE:-0} == 1 ]]; then
        echo "[$(date)] FORCE=1, removing incomplete $OUT"; rm -rf "$OUT"
    else
        echo "FATAL: $OUT exists but is incomplete (a failed or running attempt)."
        echo "       Inspect it, then re-run this task with FORCE=1 to discard and restart."
        exit 1
    fi
fi

# ---- disk headroom ----------------------------------------------------------
AVAIL_GB=$(df -BG --output=avail "$OUTROOT" | tail -1 | tr -dc '0-9')
NEED_GB=$([[ $CREATE_BAM == true ]] && echo 150 || echo 60)
[[ $AVAIL_GB -ge $NEED_GB ]] || { echo "FATAL: only ${AVAIL_GB}G free under $OUTROOT, need ~${NEED_GB}G"; exit 1; }
echo "disk free  : ${AVAIL_GB}G (need ~${NEED_GB}G)"

# keep martian's scratch on /scratch rather than the node's /tmp
export TMPDIR=$OUTROOT/.tmp_${SAMPLE}
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

# ---- count ------------------------------------------------------------------
CORES=${SLURM_CPUS_PER_TASK:-16}
MEMGB=$(( ${SLURM_MEM_PER_NODE:-131072} / 1024 - 8 ))   # leave headroom under the cgroup
echo "[$(date)] running cellranger count -> $OUT  (${CORES} cores, ${MEMGB} GB)"
cd "$OUTROOT"
"$CELLRANGER" count \
    --id="$SAMPLE" \
    --transcriptome="$REF" \
    --fastqs="$FASTQS" \
    --sample="$SAMPLE" \
    --create-bam="$CREATE_BAM" \
    --include-introns="$INCLUDE_INTRONS" \
    --localcores="$CORES" \
    --localmem="$MEMGB" \
    --disable-cell-annotation \
    --disable-ui \
    --output-dir="$OUT"

# ---- verify -----------------------------------------------------------------
echo "[$(date)] verifying outputs"
for f in outs/metrics_summary.csv outs/web_summary.html \
         outs/filtered_feature_bc_matrix/features.tsv.gz \
         outs/raw_feature_bc_matrix/features.tsv.gz; do
    [[ -s $OUT/$f ]] || { echo "FAIL: missing $OUT/$f"; exit 1; }
    echo "  OK  $f"
done

# the four added contigs must be carried through as features, or the transgene
# calls downstream are silently impossible
for g in EGFP tdTomato EBFP2 AAV_shared; do
    if zgrep -qw "$g" "$OUT/outs/raw_feature_bc_matrix/features.tsv.gz"; then
        echo "  OK  $g present in feature list"
    else
        echo "  FAIL $g MISSING from feature list"; exit 1
    fi
done

echo "--- metrics ---"
cat "$OUT/outs/metrics_summary.csv"
echo "[$(date)] done. $SAMPLE at $OUT"
du -sh "$OUT"
