#!/usr/bin/env python3
"""
Step 2.11: projection calls from the v2 matrix plus BOTH rescue tracks.

Three states are reported side by side so each track's contribution is visible
and either can be dropped:

  base   the v2 filtered matrix alone (step 2.1)
  +resc  plus BLAST-rescued molecules from the unmapped pool (steps 2.7 / 2.8)
  +anti  plus antisense molecules already aligned to the cassettes (step 2.10)

Both added tracks passed the same two controls, which is why they are folded in
at all:
  - failed-reporter: tdTomato molecules stay near zero in Ech2/6/10/11/12, where
    tdTomato failed biologically (rescue 4 vs 573; antisense 2 vs 288)
  - cellularity: rescued molecules are 72.6x enriched on cells already positive
    for the same construct; antisense reads reach called cells at 1.26x the rate
    of the SENSE reads on the same contigs

Classification is the rule from 2_1_projection_validation.R: exactly one ORF at
>=1 UMI gives that ORF's region for this mouse, two or more give `multi_ORF`, no
ORF but shared >=1 gives `shared_only`, nothing gives `unlabelled`. The
fluorophore -> region map is randomised per mouse and read from
Data/sample_metadata.csv.

Both tracks bypass a filter the mouse genes passed - BLAST-rescued reads failed
STAR's length/score cutoff, antisense reads are not counted by design. These
counts are valid for presence/absence projection labelling ONLY, never for
expression-level comparison, and the base column is kept so that boundary stays
visible.
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
ANTI = os.path.join(PROJ, "Results", "2_10_antisense_harvest")
OUTDIR = os.path.join(PROJ, "Results", "2_11_combined_projection")

ORFS = ["EGFP", "tdTomato", "EBFP2"]
REGIONS = ["BLA", "NAc", "LHr"]
STATES = ["base", "resc", "anti"]


def load_meta():
    with open(META, newline="") as fh:
        txt = fh.read().replace("\r", "")
    out = {}
    for r in csv.DictReader(txt.splitlines()):
        m = {}
        for reg in REGIONS:
            f = r[reg].strip()
            m["EBFP2" if f == "EBFP" else f] = reg
        if set(m) != set(ORFS):
            raise SystemExit("%s: fluorophore map is not a permutation of %s"
                             % (r["fastq_prefix"], ORFS))
        out[r["fastq_prefix"]] = (r["sex"], r["diet"], m)
    return out


def barcodes(sample):
    p = os.path.join(COUNTS, sample, "outs", "filtered_feature_bc_matrix",
                     "barcodes.tsv.gz")
    with gzip.open(p, "rt") as fh:
        return [l.strip() for l in fh]


def matrix_counts():
    per = defaultdict(dict)
    with gzip.open(PERCELL, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            per[r["sample"]][int(r["col"])] = {
                "EGFP": int(r["EGFP"]), "tdTomato": int(r["tdTomato"]),
                "EBFP2": int(r["EBFP2"]), "shared": int(r["shared"])}
    return per


def rescued(sample):
    p = os.path.join(RESCUE, sample + "_rescued_reads.tsv")
    seen = defaultdict(lambda: defaultdict(set))
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["in_called_cell"] == "1" and r["CB"] and r["UB"]:
                seen[r["CB"]][r["construct"].replace("_cassette", "")].add(r["UB"])
    return {b: {k: len(v) for k, v in d.items()} for b, d in seen.items()}


def antisense(sample):
    """Rows are already restricted to called cells with CB and UB."""
    p = os.path.join(ANTI, sample + "_antisense_reads.tsv")
    seen = defaultdict(lambda: defaultdict(set))
    with open(p) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            feat = "shared" if r["contig"] == "AAV_shared" \
                else r["contig"].replace("_cassette", "")
            seen[r["CB"]][feat].add(r["UB"])
    return {b: {k: len(v) for k, v in d.items()} for b, d in seen.items()}


def add(dst, extra):
    for bc, d in extra.items():
        t = dst.setdefault(bc, {})
        for k, v in d.items():
            t[k] = t.get(k, 0) + v


def classify(c, fmap):
    hot = [o for o in ORFS if c.get(o, 0) >= 1]
    if len(hot) == 1:
        return fmap[hot[0]]
    if len(hot) > 1:
        return "multi_ORF"
    return "shared_only" if c.get("shared", 0) >= 1 else "unlabelled"


def tally(cells, fmap, n_cells):
    t = defaultdict(int)
    for c in cells.values():
        t[classify(c, fmap)] += 1
        if any(c.get(o, 0) >= 1 for o in ORFS):
            t["orf_pos"] += 1
        if any(c.get(k, 0) >= 1 for k in ORFS + ["shared"]):
            t["tg_pos"] += 1
    t["assigned"] = sum(t[r] for r in REGIONS)
    t["unlabelled"] += n_cells - len(cells)
    return t


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    meta, mat = load_meta(), matrix_counts()
    rows, per = [], {}

    for sample in sorted(meta, key=lambda s: int("".join(c for c in s[3:5] if c.isdigit()))):
        sex, diet, fmap = meta[sample]
        bcs = barcodes(sample)
        n_cells = len(bcs)

        state = {"base": {bcs[col - 1]: dict(c) for col, c in mat.get(sample, {}).items()}}
        state["resc"] = {b: dict(c) for b, c in state["base"].items()}
        add(state["resc"], rescued(sample))
        state["anti"] = {b: dict(c) for b, c in state["resc"].items()}
        add(state["anti"], antisense(sample))

        per[sample] = {s: tally(state[s], fmap, n_cells) for s in STATES}
        rows.append((sample, sex, diet, n_cells))

    keys = ["tg_pos", "orf_pos"] + REGIONS + ["assigned", "multi_ORF", "shared_only"]
    with open(os.path.join(OUTDIR, "projection_three_state.tsv"), "w") as fh:
        fh.write("sample\tsex\tdiet\ttotal_cells\t"
                 + "\t".join("%s_%s" % (k, s) for k in keys for s in STATES) + "\n")
        for sample, sex, diet, n in rows:
            fh.write("%s\t%s\t%s\t%d\t" % (sample, sex, diet, n)
                     + "\t".join(str(per[sample][s][k]) for k in keys for s in STATES) + "\n")

    def tot(k, s):
        return sum(per[x[0]][s][k] for x in rows)

    print("Transgene-positive cells per animal  (base -> +BLAST rescue -> +antisense)")
    print("=" * 84)
    print("%-11s %3s %4s %8s | %26s | %26s"
          % ("sample", "sex", "diet", "cells", "transgene-positive", "ORF-positive"))
    for sample, sex, diet, n in rows:
        p = per[sample]
        print("%-11s %3s %4s %8d | %7d ->%7d ->%7d | %7d ->%7d ->%7d"
              % (sample, sex, diet, n,
                 p["base"]["tg_pos"], p["resc"]["tg_pos"], p["anti"]["tg_pos"],
                 p["base"]["orf_pos"], p["resc"]["orf_pos"], p["anti"]["orf_pos"]))
    print("%-11s %3s %4s %8d | %7d ->%7d ->%7d | %7d ->%7d ->%7d"
          % ("TOTAL", "", "", tot("tg_pos", "base") and sum(x[3] for x in rows),
             tot("tg_pos", "base"), tot("tg_pos", "resc"), tot("tg_pos", "anti"),
             tot("orf_pos", "base"), tot("orf_pos", "resc"), tot("orf_pos", "anti")))

    print()
    print("Projection assignment  (base -> +rescue -> +antisense)")
    print("=" * 84)
    print("%-11s %22s %22s %22s" % ("sample", "BLA", "NAc", "LHr"))
    for sample, sex, diet, n in rows:
        p = per[sample]
        print("%-11s %6d ->%6d ->%6d %6d ->%6d ->%6d %6d ->%6d ->%6d"
              % (sample,
                 p["base"]["BLA"], p["resc"]["BLA"], p["anti"]["BLA"],
                 p["base"]["NAc"], p["resc"]["NAc"], p["anti"]["NAc"],
                 p["base"]["LHr"], p["resc"]["LHr"], p["anti"]["LHr"]))
    print("%-11s %6d ->%6d ->%6d %6d ->%6d ->%6d %6d ->%6d ->%6d"
          % ("TOTAL",
             tot("BLA", "base"), tot("BLA", "resc"), tot("BLA", "anti"),
             tot("NAc", "base"), tot("NAc", "resc"), tot("NAc", "anti"),
             tot("LHr", "base"), tot("LHr", "resc"), tot("LHr", "anti")))

    n_all = sum(x[3] for x in rows)
    print()
    for k in ("assigned", "multi_ORF", "shared_only"):
        print("%-12s %6d -> %6d -> %6d" % (k, tot(k, "base"), tot(k, "resc"), tot(k, "anti")))
    print()
    print("assigned: %.2f%% -> %.2f%% -> %.2f%% of %d nuclei"
          % (100 * tot("assigned", "base") / n_all,
             100 * tot("assigned", "resc") / n_all,
             100 * tot("assigned", "anti") / n_all, n_all))
    print("written to", OUTDIR)


if __name__ == "__main__":
    main()
