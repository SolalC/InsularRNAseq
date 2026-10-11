#!/bin/bash
#SBATCH --job-name=hfd_backup
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=2-00:00:00
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/backup_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/backup_%j.err

# Archive the whole project to /QRISdata/Q9732 by running the project's
# own push_to_Q.sh under the scheduler.
#
# push_to_Q.sh is called, not reimplemented, so there is one definition of what
# gets copied and where. It carries no --exclude and that is deliberate
# (confirmed by the user 2026-08-29): Data/ is meant to go to the archive along
# with everything else. The payload is ~1.6 TB, of which the 767 GB of fastq and
# the two reference builds were already pushed on 2026-08-26, so this run is
# mostly incremental - but rsync still has to stat every file, hence the 2-day
# walltime rather than a few hours.
#
# This is I/O bound end to end. 2 CPU and 8 GB is enough; asking for more would
# only make it queue longer.

set -euo pipefail

PROJ=/scratch/user/uqschauq/Insular
PUSH=$PROJ/push_to_Q.sh
DEST=/QRISdata/Q9732/scratch

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"

# ---- preflight --------------------------------------------------------------
# /QRISdata is an RDM mount and is not guaranteed to be present on every compute
# node. Check it here so a missing mount fails in seconds with a readable
# message instead of rsync creating a real directory over the empty mountpoint
# and quietly filling the node's local disk.
[[ -s $PUSH ]] || { echo "FATAL: missing $PUSH"; exit 1; }
[[ -d $DEST ]] || { echo "FATAL: $DEST not visible from $(hostname) - the QRISdata mount is absent on this node"; exit 1; }
[[ -w $DEST ]] || { echo "FATAL: $DEST is not writable"; exit 1; }
mountpoint -q /QRISdata || echo "NOTE: /QRISdata is not reported as a mountpoint here; continuing because $DEST is present and writable"

echo "[$(date)] source: $(du -sh $PROJ 2>/dev/null | cut -f1) at $PROJ"
echo "[$(date)] dest before: $(du -sh $DEST 2>/dev/null | cut -f1) at $DEST"
echo "--- push_to_Q.sh ---"; cat "$PUSH"; echo "--------------------"

# ---- run --------------------------------------------------------------------
bash "$PUSH"

echo "[$(date)] dest after: $(du -sh $DEST 2>/dev/null | cut -f1)"

# ---- verify -----------------------------------------------------------------
# A second rsync pass in dry-run: anything it still wants to transfer is
# something the first pass did not land.
echo "--- residual differences (dry run) ---"
rsync -rlt --dry-run --itemize-changes "$PROJ/" "$DEST/" | head -50
echo "[$(date)] done"
