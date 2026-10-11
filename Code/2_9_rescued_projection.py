#!/usr/bin/env python3
"""
Step 2.9: fold the 2.8 rescued molecules into the per-cell transgene counts and
re-run the projection classification, so the NET effect on assigned cells is
visible rather than the gross "+780 (cell, construct) positives".

The gross figure overstates the gain. A rescued molecule landing on a cell that
is already positive for a DIFFERENT construct turns a clean single-ORF call into
`multi_ORF`, which is ambiguous and therefore unassigned. Those cells are lost
from a region, so the net change can be smaller than the gross - or negative for
an individual region.

Classification is exactly the rule 2_1_projection_validation.R uses:
  exactly one ORF >= 1 UMI -> that ORF's region for this mouse
  two or three ORFs        -> multi_ORF (ambiguous, no region)
  no ORF but shared >= 1   -> shared_only  ("transgene-positive, unassigned")
  nothing                  -> unlabelled
The fluorophore -> region map is randomised per mouse and is read per sample
from Data/sample_metadata.csv; never assume a fixed correspondence.

Rescued molecules only ever add to the three ORF counts: 2.8 kept only HSPs on
construct-diagnostic sequence, so nothing lands on the shared element.
"""

import csv
import gzip
import os
from collections import defaultdict

PROJ = "/scratch/user/uqschauq/Insular"
META = os.path.join(PROJ, "Data", "sample_metadata.csv")
COUNTS = os.path.join(PROJ, "Data", "Counts_v2")
PERCELL = os.path.join(PROJ, "Results", "2_1_projection_validation_v2",
                       "per_cell_calls_ge1.tsv.gz")
RESCUE = os.path.join(PROJ, "Results", "2_8_rescue_assign")
OUTDIR = os.path.join(PROJ, "Results", "2_9_rescued_projection")

ORFS = ["EGFP", "tdTomato", "EBFP2"]
REGIONS = ["BLA", "NAc", "LHr"]


def load_meta():
    """sample -> (sex, diet, {fluorophore: region})."""
    with open(META, newline="") as fh:
        txt = fh.read().replace("\r", "")
    rows = list(csv.DictReader(txt.splitlines()))
    out = {}
    for r in rows:
        m = {}
        for reg in REGIONS:
            f = r[reg].strip()
            m["EBFP2" if f == "EBFP" else f] = reg   # metadata says EBFP
        if set(m) != set(ORFS):
            raise SystemExit("%s: fluorophore map is not a permutation of %s: %s"
                             % (r["fastq_prefix"], ORFS, sorted(m)))
        out[r["fastq_prefix"]] = (r["sex"], r["diet"], m)
    return out


def barcodes(sample):
    p = os.path.join(COUNTS, sample, "outs", "filtered_feature_bc_matrix",
                     "barcodes.tsv.gz")
    with gzip.open(p, "rt") as fh:
        return [l.strip() for l in fh]


def matrix_counts():
    """sample -> col -> {feature: umis} from the v2 filtered matrix."""
    per = defaultdict(dict)
    with gzip.open(PERCELL, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            per[r["sample"]][int(r["col"])] = {
                "EGFP": int(r["EGFP"]), "tdTomato": int(r["tdTomato"]),
                "EBFP2": int(r["EBFP2"]), "shared": int(r["shared"])}
    return per


def rescued_molecules(sample):
    """barcode -> {ORF: n distinct UMIs}, called cells only."""
    p = os.path.join(RESCUE, sample + "_rescued_reads.tsv")
    seen = defaultdict(lambda: defaultdict(set))
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["in_called_cell"] != "1" or not r["CB"] or not r["UB"]:
                continue
            seen[r["CB"]][r["construct"].replace("_cassette", "")].add(r["UB"])
    return {bc: {o: len(u) for o, u in d.items()} for bc, d in seen.items()}


def classify(counts, fmap):
    """counts = {feature: umis} -> class label."""
    hot = [o for o in ORFS if counts.get(o, 0) >= 1]
    if len(hot) == 1:
        return fmap[hot[0]]
    if len(hot) > 1:
        return "multi_ORF"
    if counts.get("shared", 0) >= 1:
        return "shared_only"
    return "unlabelled"


def tally(cells, fmap):
    t = defaultdict(int)
    for c in cells.values():
        t[classify(c, fmap)] += 1
        if any(c.get(o, 0) >= 1 for o in ORFS):
            t["orf_pos"] += 1
        if any(c.get(k, 0) >= 1 for k in ORFS + ["shared"]):
            t["transgene_pos"] += 1
    t["assigned"] = sum(t[r] for r in REGIONS)
    return t


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    meta = load_meta()
    mat = matrix_counts()

    hdr = ("sample sex diet total_cells "
           "tg_pos_before tg_pos_after orf_pos_before orf_pos_after "
           "BLA_before BLA_after NAc_before NAc_after LHr_before LHr_after "
           "assigned_before assigned_after multi_before multi_after "
           "shared_only_before shared_only_after").split()
    rows = []

    for sample in sorted(meta, key=lambda s: int("".join(ch for ch in s[3:5] if ch.isdigit()))):
        sex, diet, fmap = meta[sample]
        bcs = barcodes(sample)
        n_cells = len(bcs)

        before = {}
        for col, c in mat.get(sample, {}).items():
            before[bcs[col - 1]] = dict(c)

        after = {bc: dict(c) for bc, c in before.items()}
        for bc, add in rescued_molecules(sample).items():
            d = after.setdefault(bc, {})
            for o, n in add.items():
                d[o] = d.get(o, 0) + n

        # Cells absent from the tables have all-zero counts; they can only be
        # "unlabelled", so they need no entry - but they do count in the total.
        tb, ta = tally(before, fmap), tally(after, fmap)
        zero = n_cells - len(before)
        tb["unlabelled"] += zero
        ta["unlabelled"] += n_cells - len(after)

        rows.append([sample, sex, diet, n_cells,
                     tb["transgene_pos"], ta["transgene_pos"],
                     tb["orf_pos"], ta["orf_pos"],
                     tb["BLA"], ta["BLA"], tb["NAc"], ta["NAc"],
                     tb["LHr"], ta["LHr"],
                     tb["assigned"], ta["assigned"],
                     tb["multi_ORF"], ta["multi_ORF"],
                     tb["shared_only"], ta["shared_only"]])

    with open(os.path.join(OUTDIR, "projection_before_after.tsv"), "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    def col(i):
        return sum(r[i] for r in rows)

    print("Transgene-positive cells per animal (>=1 UMI on any transgene feature)")
    print("=" * 78)
    print("%-11s %3s %4s %8s | %9s %9s %6s | %8s %8s"
          % ("sample", "sex", "diet", "cells", "tg_before", "tg_after", "delta",
             "orf_bef", "orf_aft"))
    for r in rows:
        print("%-11s %3s %4s %8d | %9d %9d %+6d | %8d %8d"
              % (r[0], r[1], r[2], r[3], r[4], r[5], r[5] - r[4], r[6], r[7]))
    print("%-11s %3s %4s %8d | %9d %9d %+6d | %8d %8d"
          % ("TOTAL", "", "", col(3), col(4), col(5), col(5) - col(4), col(6), col(7)))

    print()
    print("Projection assignment, before -> after rescue")
    print("=" * 78)
    print("%-11s %13s %13s %13s | %13s %9s %11s"
          % ("sample", "BLA", "NAc", "LHr", "assigned", "multi", "shared_only"))
    for r in rows:
        print("%-11s %6d ->%5d %6d ->%5d %6d ->%5d | %6d ->%5d %4d ->%4d %5d ->%5d"
              % (r[0], r[8], r[9], r[10], r[11], r[12], r[13],
                 r[14], r[15], r[16], r[17], r[18], r[19]))
    print("%-11s %6d ->%5d %6d ->%5d %6d ->%5d | %6d ->%5d %4d ->%4d %5d ->%5d"
          % ("TOTAL", col(8), col(9), col(10), col(11), col(12), col(13),
             col(14), col(15), col(16), col(17), col(18), col(19)))
    print()
    print("net change in assigned cells: %+d (%.2f%% -> %.2f%% of nuclei)"
          % (col(15) - col(14), 100 * col(14) / col(3), 100 * col(15) / col(3)))
    print("written to", OUTDIR)


if __name__ == "__main__":
    main()
