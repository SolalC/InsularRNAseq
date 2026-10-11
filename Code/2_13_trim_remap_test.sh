#!/bin/bash
#SBATCH --job-name=hfd_trimremap
#SBATCH --partition=general
#SBATCH --array=1-12
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/trimremap_%A_%a.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/trimremap_%A_%a.err

# Step 2.13: WOULD BETTER READ QUALITY ACTUALLY HELP THE MAPPING?
#
# FastQC (step 0) describes the reads; it cannot say whether fixing what it
# describes would recover anything. This measures it directly.
#
# Take the reads cellranger failed to map, trim them with fastp (poly-G, poly-X,
# adapters, quality), and re-align them to the SAME v2 reference. Every read that
# now aligns is a read that better quality handling would have rescued; every
# read that still fails is one that trimming cannot help.
#
# The instrument is a NextSeq 2000 (VH00807), two-colour chemistry, where G means
# "no signal" - so dark cycles produce poly-G tails. That is the artefact most
# likely to be costing alignments here, since the unmapped records carry
# `uT:A:1` (too short) alongside otherwise perfect partial alignments, e.g.
# `AS:i:51 nM:i:0` on a 90 bp read.
#
# STAR is run at DEFAULTS (outFilterScoreMinOverLread / MatchNminOverLread 0.66).
# That is not an arbitrary choice: 0.66 * 90 = 59.4, and the rejected records show
# alignments at 51 bp surviving to `uT:A:1`, so the default reproduces the bar
# cellranger actually applied. The comparison is therefore like for like, and
# because the filters are RELATIVE to read length, a trimmed read needs
# proportionally less matched sequence - which is precisely how trimming rescues
# an alignment.

set -euo pipefail

module purge
# fastp 0.23.2 is built against gcc/11.3.0 and samtools/1.18 against gcc/12.3.0;
# this site forbids same-name module autoswap, so both must come from the same
# toolchain. samtools 1.16.1 is the gcc-11.3.0 build and does everything needed
# here (indexed region query + view).
module load samtools/1.16.1-gcc-11.3.0
module load fastp/0.23.2-gcc-11.3.0

PROJ=/scratch/user/uqschauq/Insular
STAR=$PROJ/Soft/cellranger-10.1.0/lib/bin/STAR
IDX=$PROJ/Data/Reference/GRCm39_cassette_v2
MANIFEST=$PROJ/Data/sample_metadata.csv
OUT=$PROJ/Results/2_13_trim_remap
CORES=${SLURM_CPUS_PER_TASK:-8}

mkdir -p "$OUT"
TMP=${TMPDIR:-/tmp}/trimremap.$SLURM_ARRAY_JOB_ID.$SLURM_ARRAY_TASK_ID
mkdir -p "$TMP"
trap 'rm -rf "$TMP"' EXIT

# sample_metadata.csv is CRLF; strip \r and look the column up by header name.
SAMPLE=$(awk -F',' -v n="$SLURM_ARRAY_TASK_ID" '
  NR==1 { for (i=1;i<=NF;i++) { h=$i; gsub(/\r/,"",h); if (h=="fastq_prefix") c=i }; next }
  NR==n+1 { v=$c; gsub(/\r/,"",v); print v }' "$MANIFEST")
[[ -n $SAMPLE ]] || { echo "FATAL: no sample for array index $SLURM_ARRAY_TASK_ID"; exit 1; }

BAM=$PROJ/Data/Counts_v2/$SAMPLE/outs/possorted_genome_bam.bam
[[ -s $BAM ]]        || { echo "FATAL: no BAM at $BAM"; exit 1; }
[[ -x $STAR ]]       || { echo "FATAL: no STAR at $STAR"; exit 1; }
[[ -s $IDX/star/SA ]] || { echo "FATAL: no v2 STAR index under $IDX"; exit 1; }

echo "[$(date)] $SAMPLE: extracting unmapped reads"
# Column 11 is the quality string; unmapped records keep the original sequence
# and qualities, so this reconstructs the read as it came off the instrument.
samtools view -@ 2 "$BAM" '*' \
  | awk '{ print "@" $1 "\n" $10 "\n+\n" $11 }' > "$TMP/unmapped.fq"
IN_READS=$(( $(wc -l < "$TMP/unmapped.fq") / 4 ))
echo "[$(date)] $SAMPLE: $IN_READS unmapped reads"

echo "[$(date)] $SAMPLE: fastp"
# --trim_poly_g is the NextSeq two-colour fix; --trim_poly_x catches polyA/T
# runoff. Length floor of 25 keeps anything STAR could conceivably place.
fastp \
  --in1 "$TMP/unmapped.fq" --out1 "$TMP/trimmed.fq" \
  --trim_poly_g --trim_poly_x \
  --cut_tail --cut_tail_mean_quality 20 \
  --length_required 25 \
  --thread "$CORES" \
  --json "$OUT/${SAMPLE}_fastp.json" \
  --html "$OUT/${SAMPLE}_fastp.html" \
  2> "$OUT/${SAMPLE}_fastp.log"
OUT_READS=$(( $(wc -l < "$TMP/trimmed.fq") / 4 ))
echo "[$(date)] $SAMPLE: $OUT_READS reads survive trimming"

echo "[$(date)] $SAMPLE: re-aligning with STAR at defaults"
"$STAR" \
  --runThreadN "$CORES" \
  --genomeDir "$IDX/star" \
  --readFilesIn "$TMP/trimmed.fq" \
  --outFileNamePrefix "$TMP/remap_" \
  --outSAMtype SAM \
  --outStd SAM \
  --outSAMunmapped None \
  | awk -v s="$SAMPLE" -v inr="$IN_READS" -v outr="$OUT_READS" '
      /^@/ { next }
      { if (and($2,4)) next
        n++
        c=$3
        if (c=="EGFP_cassette" || c=="tdTomato_cassette" || c=="EBFP2_cassette" || c=="AAV_shared") tg[c]++
      }
      END {
        printf "metric\tvalue\n"                      > OUTF
        printf "sample\t%s\n", s                      > OUTF
        printf "unmapped_reads_in\t%d\n", inr         > OUTF
        printf "reads_after_fastp\t%d\n", outr        > OUTF
        printf "newly_mapped\t%d\n", n                > OUTF
        printf "newly_mapped_pct_of_unmapped\t%.4f\n", (inr ? 100*n/inr : 0) > OUTF
        for (c in tg) printf "newly_mapped_%s\t%d\n", c, tg[c] > OUTF
      }' OUTF="$OUT/${SAMPLE}_remap_summary.tsv"

cp "$TMP/remap_Log.final.out" "$OUT/${SAMPLE}_STAR_Log.final.out" 2>/dev/null || true
echo "[$(date)] $SAMPLE: result"
cat "$OUT/${SAMPLE}_remap_summary.tsv"
echo "[$(date)] done"
