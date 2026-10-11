#!/bin/bash
#SBATCH --job-name=hfd_mkref_v2
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/mkrefv2_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/mkrefv2_%j.err

# Build the HFDseq v2 reference: 10x prebuilt mouse GRCm39-2024-A plus four
# cassette-derived transgene contigs.
#
# v1 held three ORF fragments and one 799 bp "shared" 3' element cut from p24.
# That element is not shared - p131 matches 395/799 bp, p140 589/799 - so the 3'
# end of the tdTomato and EBFP2 constructs was absent from the reference and
# 7,952 diagnostic reads (+32%) could not be counted. See
# Results/2_4_cassette_report/FINDINGS.md.
#
# v2 contigs are generated here by 1_mkref_cassette_v2.py from the reconstructed
# ITR-to-ITR cassettes, so every distinct sequence is present exactly once:
# <reporter>_cassette carries that construct's diagnostic bases (the rest masked
# to N) and AAV_shared carries the blocks common to more than one construct.

set -euo pipefail

PROJ=/scratch/user/uqschauq/Insular
CELLRANGER=$PROJ/Soft/cellranger-10.1.0/bin/cellranger
REFDIR=$PROJ/Data/Reference
BASE=$REFDIR/refdata-gex-GRCm39-2024-A
TG=$REFDIR/transgenes_cassette_v2
BUILD=$REFDIR/build
OUT=$REFDIR/GRCm39_cassette_v2

mkdir -p "$BUILD" "$PROJ/Code/logs"

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"

# ---- generate the v2 transgene contigs from the cassettes -------------------
echo "[$(date)] building v2 transgene contigs"
python3 "$PROJ/Code/1_mkref_cassette_v2.py"

# ---- checks on inputs -------------------------------------------------------
[[ -x $CELLRANGER ]] || { echo "FATAL: cellranger not executable at $CELLRANGER"; exit 1; }
[[ -d $BASE ]]       || { echo "FATAL: 10x reference not found at $BASE"; exit 1; }

GENOME_FA=$BASE/fasta/genome.fa
[[ -s $GENOME_FA ]]  || { echo "FATAL: missing $GENOME_FA"; exit 1; }

# 2024-A ships genes.gtf gzipped; tolerate either form
if   [[ -s $BASE/genes/genes.gtf.gz ]]; then GENES_SRC=$BASE/genes/genes.gtf.gz; GZ=yes
elif [[ -s $BASE/genes/genes.gtf    ]]; then GENES_SRC=$BASE/genes/genes.gtf;    GZ=no
else echo "FATAL: no genes.gtf(.gz) under $BASE/genes/"; exit 1; fi
echo "base fasta : $GENOME_FA"
echo "base gtf   : $GENES_SRC"

for f in transgenes.fa transgenes.gtf; do
    [[ -s $TG/$f ]] || { echo "FATAL: missing $TG/$f"; exit 1; }
done

# transgene contig names must appear in both files, and must not collide with mouse contigs
comm -3 <(grep '^>' "$TG/transgenes.fa" | sed 's/^>//' | sort) \
        <(cut -f1 "$TG/transgenes.gtf" | sort -u) > "$BUILD/name_mismatch.txt"
if [[ -s $BUILD/name_mismatch.txt ]]; then
    echo "FATAL: transgene contig names differ between FASTA and GTF:"; cat "$BUILD/name_mismatch.txt"; exit 1
fi

grep '^>' "$GENOME_FA" | sed 's/^>//' | cut -d' ' -f1 | sort > "$BUILD/mouse_contigs.txt"
if comm -12 "$BUILD/mouse_contigs.txt" \
            <(grep '^>' "$TG/transgenes.fa" | sed 's/^>//' | sort) | grep -q .; then
    echo "FATAL: a transgene contig name collides with a mouse contig name"; exit 1
fi
echo "input checks passed"

# ---- assemble combined fasta + gtf -----------------------------------------
COMB_FA=$BUILD/genome_plus_transgenes.fa
COMB_GTF=$BUILD/genes_plus_transgenes.gtf

echo "[$(date)] assembling combined FASTA"
cat "$GENOME_FA" > "$COMB_FA"
# guarantee a newline between the two files before appending
if [[ $(tail -c1 "$COMB_FA" | wc -l) -eq 0 ]]; then echo >> "$COMB_FA"; fi
cat "$TG/transgenes.fa" >> "$COMB_FA"

echo "[$(date)] assembling combined GTF"
if [[ $GZ == yes ]]; then zcat "$GENES_SRC" > "$COMB_GTF"; else cp "$GENES_SRC" "$COMB_GTF"; fi
cat "$TG/transgenes.gtf" >> "$COMB_GTF"

echo "combined FASTA contigs: $(grep -c '^>' "$COMB_FA")  (base $(grep -c '^>' "$GENOME_FA") + $(grep -c '^>' "$TG/transgenes.fa"))"
echo "combined GTF lines    : $(wc -l < "$COMB_GTF")"

# ---- build ------------------------------------------------------------------
rm -rf "$OUT"
echo "[$(date)] running cellranger mkref -> $OUT"
cd "$BUILD"
"$CELLRANGER" mkref \
    --genome=GRCm39_cassette_v2 \
    --fasta="$COMB_FA" \
    --genes="$COMB_GTF" \
    --nthreads="${SLURM_CPUS_PER_TASK:-8}" \
    --memgb=60 \
    --ref-version=GRCm39-2024-A-cassette-v2 \
    --output-dir="$OUT"

# ---- confirm the transgenes survived into the built reference ---------------
echo "[$(date)] verifying built reference"
BUILT_GTF=$OUT/genes/genes.gtf.gz
[[ -s $BUILT_GTF ]] || BUILT_GTF=$OUT/genes/genes.gtf
for g in EGFP tdTomato EBFP2 AAV_shared; do
    if zgrep -q "gene_id \"$g\"" "$BUILT_GTF" 2>/dev/null || grep -q "gene_id \"$g\"" "$BUILT_GTF"; then
        echo "  OK  $g present in built annotation"
    else
        echo "  FAIL $g MISSING from built annotation"; exit 1
    fi
done

echo "[$(date)] done. Reference at $OUT"
du -sh "$OUT"
