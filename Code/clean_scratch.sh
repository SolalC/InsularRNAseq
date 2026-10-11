#!/bin/bash
#SBATCH --job-name=clean_scratch
#SBATCH --account=a_imb_cpdg
#SBATCH --partition=general
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=8G --time=06:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/clean_scratch_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/clean_scratch_%j.err

# Housekeeping: free the scratch space once push_to_Q.sh has archived the
# project to /QRISdata/Q9732/scratch. Keeps only Code/ and
# Data/Seurat/HFDseq_v2_annotated.rds; everything else under the project root
# is deleted. Nothing under /QRISdata is written.
#
# Guard: before deleting anything, the push's own rsync command is re-run as a
# dry run. If it would transfer anything at all, the archive is not complete
# and the script aborts. Code/ is excluded from that check because it is kept
# (and this script plus its logs are newer than the last push).

set -euo pipefail
PROJ=/scratch/user/uqschauq/Insular
ARCHIVE=/QRISdata/Q9732/scratch
KEEP_RDS=HFDseq_v2_annotated.rds

ts() { date '+%F %T'; }

REAL=$(realpath "$PROJ")
[[ "$REAL" == /QRISdata* ]] && { echo "ABORT: $PROJ resolves to $REAL" >&2; exit 1; }
[[ -f "$PROJ/Data/Seurat/$KEEP_RDS" ]] || { echo "ABORT: $KEEP_RDS not found" >&2; exit 1; }

echo "$(ts) dry-run rsync $PROJ/ -> $ARCHIVE/"
pending=$(rsync -rlt --dry-run --itemize-changes --exclude=/Code/ "$PROJ/" "$ARCHIVE/")
if [[ -n "$pending" ]]; then
    echo "ABORT: archive is not up to date; rsync would still transfer:" >&2
    echo "$pending" | head -100 >&2
    exit 1
fi
echo "$(ts) archive complete, nothing pending"

shopt -s dotglob nullglob

delete() {
    echo "$(ts) rm $1"
    rm -rf --one-file-system -- "$1"
}

for p in "$PROJ"/*; do
    case "${p##*/}" in Code|Data) continue ;; esac
    delete "$p"
done
for p in "$PROJ"/Data/*; do
    [[ "${p##*/}" == Seurat ]] && continue
    delete "$p"
done
for p in "$PROJ"/Data/Seurat/*; do
    [[ "${p##*/}" == "$KEEP_RDS" ]] && continue
    delete "$p"
done

echo "$(ts) done; remaining:"
find "$PROJ" -maxdepth 3 -not -path "$PROJ/Code/*"
ls -l "$PROJ/Data/Seurat/$KEEP_RDS"
du -sh "$PROJ"
