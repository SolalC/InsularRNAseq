#!/bin/bash
#SBATCH --job-name=hfd_aggr
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/aggr_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/aggr_%j.err

# Step 3: cellranger aggr - merge the 12 per-mouse v2 count matrices into one.
#
# Input is Data/Counts_v2 (the cassette_v2 reference), NOT Data/Counts (v1).
# Transgene features are therefore EGFP_cassette / tdTomato_cassette /
# EBFP2_cassette / AAV_shared.
#
# --normalize=none (user decision 2026-08-26). The default, 'mapped', subsamples
# reads to equalise median depth across libraries. Projection calls here are
# presence/absence at a 1-UMI threshold, so subsampling can delete the single
# transgene UMI a call rests on. Depth normalisation happens downstream in Seurat
# instead, where it can be applied to the mouse genes and withheld from the
# transgene features.
#
# All 12 mice are included. Ech6MHFD is QC-normal (30,125 nuclei, the highest of
# the 12) and failed only in the transgene arm (10 of 30,125 transgene-positive);
# that exclusion belongs to the projection analysis in Seurat, not to the
# transcriptome matrix.
#
# The aggregation CSV is generated here from Data/sample_metadata.csv so there is
# one source of truth for sample names and covariates, the same rule 2_count*.sh
# follows.

set -euo pipefail

PROJ=/scratch/user/uqschauq/Insular
CELLRANGER=$PROJ/Soft/cellranger-10.1.0/bin/cellranger
COUNTROOT=$PROJ/Data/Counts_v2
MANIFEST=$PROJ/Data/sample_metadata.csv
OUTROOT=$PROJ/Data/Aggr_v2
RUNID=HFDseq_aggr_v2
CSV=$OUTROOT/aggregation_input.csv

# Secondary analysis (PCA / graph clustering / t-SNE / UMAP over ~324k nuclei) is
# the expensive part of aggr and nothing downstream reads it: clustering is done
# in Seurat. Set to false to let cellranger run it anyway.
NOSECONDARY=true

mkdir -p "$OUTROOT" "$PROJ/Code/logs"

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"

# ---- preflight --------------------------------------------------------------
[[ -x $CELLRANGER ]] || { echo "FATAL: cellranger not executable at $CELLRANGER"; exit 1; }
[[ -s $MANIFEST ]]   || { echo "FATAL: missing manifest $MANIFEST"; exit 1; }

# ---- build the aggregation CSV ----------------------------------------------
# sample_metadata.csv has CRLF line endings and column order is not guaranteed,
# so strip \r and look every column up by header name.
#
# Extra columns beyond sample_id/molecule_h5 are carried through to
# outs/aggregation.csv as per-library metadata, ready to join onto the Seurat
# object. One name is deliberately avoided: a column called 'batch' switches on
# cellranger's chemistry batch correction, which would alter the counts this run
# exists to keep raw.
awk -v countroot="$COUNTROOT" '
  BEGIN { FS = ","; OFS = "," }
  {
    sub(/\r$/, "")
  }
  NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    for (k in col) if (tolower(k) == "batch") {
      print "FATAL: manifest has a column named batch" > "/dev/stderr"; exit 1
    }
    split("sample_id fastq_prefix diet sex BLA NAc LHr fastq_ok", need, " ")
    for (i in need) if (!(need[i] in col)) {
      print "FATAL: manifest lacks column " need[i] > "/dev/stderr"; exit 1
    }
    print "sample_id", "molecule_h5", "mouse", "diet", "sex", "BLA", "NAc", "LHr"
    next
  }
  $(col["fastq_ok"]) == "True" {
    prefix = $(col["fastq_prefix"])
    h5 = countroot "/" prefix "/outs/molecule_info.h5"
    print prefix, h5, $(col["sample_id"]), $(col["diet"]), $(col["sex"]), \
          $(col["BLA"]), $(col["NAc"]), $(col["LHr"])
    n++
  }
  END {
    if (n != 12) { print "FATAL: " n " samples in manifest, expected 12" > "/dev/stderr"; exit 1 }
  }
' "$MANIFEST" > "$CSV"

echo "--- aggregation CSV ---"
cat "$CSV"

# Every molecule_info.h5 must exist before cellranger is started, so a missing
# input fails in seconds rather than part way through the merge.
missing=0
while IFS=, read -r sid h5 rest; do
  [[ $sid == "sample_id" ]] && continue
  [[ -s $h5 ]] || { echo "FATAL: missing $h5"; missing=1; }
done < "$CSV"
[[ $missing -eq 0 ]] || exit 1
echo "[$(date)] all 12 molecule_info.h5 present"

# ---- run --------------------------------------------------------------------
# cellranger writes its pipestance into the current directory unless told
# otherwise; --output-dir keeps it inside Data/Aggr_v2. Note the closing banner
# prints a path derived from --id, which is wrong when --output-dir is used
# (the same cosmetic mismatch 1_mkref.sh showed) - trust $OUTROOT/$RUNID.
cd "$OUTROOT"

ARGS=(--id "$RUNID"
      --csv "$CSV"
      --normalize none
      --output-dir "$OUTROOT/$RUNID"
      --localcores "${SLURM_CPUS_PER_TASK:-16}"
      --localmem 120
      --disable-ui)
[[ $NOSECONDARY == true ]] && ARGS+=(--nosecondary)

echo "[$(date)] cellranger aggr ${ARGS[*]}"
"$CELLRANGER" aggr "${ARGS[@]}"

# ---- verify -----------------------------------------------------------------
OUTS=$OUTROOT/$RUNID/outs
for f in filtered_feature_bc_matrix.h5 raw_feature_bc_matrix.h5 aggregation.csv web_summary.html; do
  [[ -s $OUTS/$f ]] || { echo "FATAL: missing $OUTS/$f"; exit 1; }
done
[[ -e $OUTROOT/$RUNID/_lock ]] && { echo "FATAL: pipestance still locked"; exit 1; }
echo "--- outs/aggregation.csv ---"
cat "$OUTS/aggregation.csv"
echo "[$(date)] done: $OUTS"
