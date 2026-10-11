#!/bin/bash
#SBATCH --job-name=inst_cellbender
#SBATCH --partition=general
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --constraint=epyc3
#SBATCH --output=/scratch/user/uqschauq/Insular/Code/logs/cellbender_install_%j.out
#SBATCH --error=/scratch/user/uqschauq/Insular/Code/logs/cellbender_install_%j.err

# Housekeeping: install CellBender into a project-local venv for step 4.
#
# The site module cellbender/0.3.0 is not usable here for two reasons: it is
# pinned to 0.3.0, which the upstream README advises against (0.3.2+ is the
# recommended line), and it lives in a read-only conda env under /sw/local so it
# cannot be upgraded in place.
#
# Version is pinned to 0.3.2 -- the version the user named. 0.4.0 exists on PyPI
# but is a newer line than the one the upstream docs recommend; bump CB_VERSION
# if that changes.
#
# torch comes from the default PyPI wheel, which bundles its own CUDA runtime,
# so no cuda module is loaded or needed. The GPU is only touched at RUN time.
#
# --constraint=epyc3 IS LOAD-BEARING. A venv hardcodes the absolute path of the
# interpreter that created it, and the module tree is per-architecture
# (/sw/auto/rocky9a/<arch>/software/Python/...). The first build of this venv
# landed on bun151, an epyc4 node, so bin/python3 pointed at the epyc4 binary;
# every gpu_cuda node except the two l40s boxes is epyc3, and the epyc4
# interpreter died there with SIGILL (bash reported 132) the instant it was
# invoked, taking all 12 array tasks down in ~3 seconds with an empty stderr.
# epyc3 is the lowest common denominator and runs forward on epyc4/epyc5, so
# build here and the venv works on every node. Same reasoning as
# Code/Makevars.portable for the R library.

set -euo pipefail

PROJ=/scratch/user/uqschauq/Insular
CB_VERSION=0.4.0
# 0.4.0, not 0.3.2. CellBender 0.3.2 declares a bare `torch` with no upper
# bound and is BROKEN against any modern torch: save_checkpoint() calls
# torch.save() on the whole pyro model object, which hits "cannot pickle
# 'weakref' object" (torch's nn.Module hook machinery holds weakrefs from 2.x
# on). Every checkpoint save then fails, and ~32 min later the run dies in
# load_or_compute_posterior_and_save() on
#   AssertionError: Checkpoint file ... does not exist
# because the posterior step hard-requires that tarball. Measured on this
# cluster: 0.3.2 fails identically with torch 2.13.0+pyro 1.9.1, with torch
# 2.3.1+pyro 1.9.1, and with torch 2.3.1+pyro 1.8.6. Going older on torch would
# fix pickling but drop H100 sm_90 support, so it is not a way out.
# TORCH_VERSION is left empty so CellBender resolves its own torch; set it only
# to pin a specific wheel.
TORCH_VERSION=""    # empty = let CellBender resolve its own torch
VENV=$PROJ/Soft/cellbender-${CB_VERSION}/venv

module purge
module load python/3.10.4-gcccore-11.3.0

# Do not let ~/.local site-packages leak into the venv.
export PYTHONNOUSERSITE=1

echo "[$(date)] host=$(hostname) job=${SLURM_JOB_ID:-none}"
echo "python: $(python3 -V 2>&1)  ($(command -v python3))"
echo "target: $VENV"

if [[ -d "$VENV" ]]; then
  echo "FATAL: $VENV already exists -- remove it deliberately before reinstalling" >&2
  exit 1
fi

mkdir -p "$(dirname "$VENV")"
python3 -m venv "$VENV"
# shellcheck disable=SC1091
source "$VENV/bin/activate"

python -m pip install --upgrade pip setuptools wheel
python -m pip install "cellbender==${CB_VERSION}"
# Installed after cellbender so the pin, if set, wins over its requirement.
if [[ -n "$TORCH_VERSION" ]]; then
  python -m pip install "torch==${TORCH_VERSION}"
fi

echo "=== verification ==="
python - <<'PY'
import importlib.metadata as md
import torch, cellbender
print("cellbender :", md.version("cellbender"))
print("torch      :", torch.__version__)
print("torch cuda :", torch.version.cuda, "| compiled with cuda:", torch.backends.cuda.is_built())
print("pyro       :", md.version("pyro-ppl"))
# get_arch_list() returns [] without a CUDA context, so it is meaningless on
# the CPU build node and is NOT checked here. sm_90 support is verified on a
# GPU node by Code/smoke_cellbender.sh instead.
print("built for cuda :", torch.version.cuda, "(sm_90 checked on GPU by the smoke test)")
PY

# --help exercises the full CLI import path; cuda is not required for this.
"$VENV/bin/cellbender" remove-background --help > /dev/null
echo "cellbender CLI OK: $("$VENV/bin/cellbender" --version 2>&1 | head -1)"

echo "[$(date)] done -- activate with: source $VENV/bin/activate"
