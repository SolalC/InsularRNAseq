#!/usr/bin/env python3
"""
Are transgene-positive nuclei being lost by cellranger's cell calling?

Reads the RAW matrix (all barcodes) on stdin as MatrixMarket triplets, and in a
single pass accumulates, per barcode: the total UMI count and the counts on the
four transgene rows. Barcodes are then split by whether cellranger called them
cells, so a population of transgene-positive barcodes sitting below the cell
threshold would show up.

Also gives the ambient-transgene estimate: transgene UMIs in clearly-empty
droplets (low total UMI), which is what a 1-UMI presence/absence call has to be
specific against.

argv: <sample> <raw_features.tsv.gz> <raw_barcodes.tsv.gz> <filtered_barcodes.tsv.gz>
"""
import sys, gzip, collections

sample, f_feat, f_bc, f_fbc = sys.argv[1:5]
ORFS = ("EGFP", "tdTomato", "EBFP2")
SHARED = "AAV_WPRE_SV40pA_shared"

# transgene row indices (1-based), looked up rather than hardcoded
rows = {}
with gzip.open(f_feat, "rt") as fh:
    for i, line in enumerate(fh, 1):
        fid = line.split("\t")[0]
        if fid in ORFS or fid == SHARED:
            rows[i] = fid
assert len(rows) == 4, f"{sample}: found {len(rows)} transgene rows, expected 4"

barcodes = [l.strip() for l in gzip.open(f_bc, "rt")]
called = set(l.strip() for l in gzip.open(f_fbc, "rt"))

total = collections.defaultdict(int)          # col -> total UMI
tg = collections.defaultdict(dict)            # col -> {gene: count}

for line in sys.stdin:
    if line[0] == "%":
        continue
    a, b, c = line.split()
    col = int(b); cnt = int(c)
    total[col] += cnt
    r = int(a)
    if r in rows:
        tg[col][rows[r]] = cnt

def orf_umi(d):    return sum(d.get(g, 0) for g in ORFS)
def shared_umi(d): return d.get(SHARED, 0)

res = collections.Counter()
res["barcodes_total"] = len(barcodes)
res["cells_called"] = len(called)

# transgene-positive barcodes, split by cell-calling status
missed = []          # non-cell barcodes with an ORF UMI, with their total UMI
ambient_orf = 0      # ORF UMIs in clearly-empty droplets (<100 total UMI)
ambient_drops = 0

for col, d in tg.items():
    bc = barcodes[col - 1]
    is_cell = bc in called
    o, s = orf_umi(d), shared_umi(d)
    tot = total[col]
    if o > 0:
        res["orf_pos_cell" if is_cell else "orf_pos_noncell"] += 1
    if s > 0:
        res["shared_pos_cell" if is_cell else "shared_pos_noncell"] += 1
    if not is_cell and o > 0:
        missed.append((tot, o, s, bc))
    if tot < 100:
        ambient_drops += 1
        ambient_orf += o

# how many excluded barcodes look like real cells rather than empty droplets
for cut in (100, 500, 1000, 2000):
    res[f"noncell_orfpos_totalUMI_ge_{cut}"] = sum(1 for t, o, s, b in missed if t >= cut)

missed.sort(reverse=True)
print(f"#sample\t{sample}")
for k in ("barcodes_total", "cells_called", "orf_pos_cell", "orf_pos_noncell",
          "shared_pos_cell", "shared_pos_noncell",
          "noncell_orfpos_totalUMI_ge_100", "noncell_orfpos_totalUMI_ge_500",
          "noncell_orfpos_totalUMI_ge_1000", "noncell_orfpos_totalUMI_ge_2000"):
    print(f"{k}\t{res[k]}")
print(f"empty_droplets_lt100UMI\t{ambient_drops}")
print(f"ambient_ORF_UMIs_in_those\t{ambient_orf}")
print(f"#top_noncell_orfpos\ttotalUMI\torfUMI\tsharedUMI\tbarcode")
for t, o, s, b in missed[:20]:
    print(f"#top\t{t}\t{o}\t{s}\t{b}")
