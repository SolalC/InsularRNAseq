#!/usr/bin/env python3
"""
Step 2.12: are the projection calls unique per nucleus?

`multi_ORF` nuclei - two or three reporters over threshold - are currently
excluded from every region, so this asks how many there are, which region pairs
they represent, and whether there are MORE of them than coincidence predicts.

The distinction matters biologically. A nucleus carrying two reporters is one of:
  - a collateralising neuron that genuinely projects to both targets (real, and
    excluding it discards a real biological class)
  - a droplet multiplet, two nuclei sharing a barcode
  - coincident ambient/index-hopping molecules on one nucleus
  - EGFP/EBFP2 cross-mapping (they are 97% identical over the ORF)

The first is the only one worth keeping, and the way to tell is the expected
count under independence. Random multiplets and independent ambient both put
reporters together at the product of their marginal rates:
    expected(A and B) = n_cells * p_A * p_B
A genuine collateral population sits ABOVE that; coincidence sits at it. The
ratio is reported per reporter pair and per region pair.

Counts are the same three states as 2.11 (base / +BLAST rescue / +antisense) so
it is visible whether the rescue tracks manufactured the overlaps.
"""

import csv
import gzip
import os
from collections import defaultdict
from itertools import combinations

PROJ = "/scratch/user/uqschauq/Insular"
META = os.path.join(PROJ, "Data", "sample_metadata.csv")
COUNTS = os.path.join(PROJ, "Data", "Counts_v2")
PERCELL = os.path.join(PROJ, "Results", "2_1_projection_validation_v2",
                       "per_cell_calls_ge1.tsv.gz")
RESCUE = os.path.join(PROJ, "Results", "2_8_rescue_assign")
ANTI = os.path.join(PROJ, "Results", "2_10_antisense_harvest")
OUTDIR = os.path.join(PROJ, "Results", "2_12_multiorf_analysis")

ORFS = ["EGFP", "tdTomato", "EBFP2"]
REGIONS = ["BLA", "NAc", "LHr"]


def load_meta():
    with open(META, newline="") as fh:
        txt = fh.read().replace("\r", "")
    out = {}
    for r in csv.DictReader(txt.splitlines()):
        m = {}
        for reg in REGIONS:
            f = r[reg].strip()
            m["EBFP2" if f == "EBFP" else f] = reg
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
            per[r["sample"]][int(r["col"])] = {o: int(r[o]) for o in ORFS}
    return per


def extra(path, sample, fname, keymap):
    seen = defaultdict(lambda: defaultdict(set))
    with open(os.path.join(path, sample + fname)) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = keymap(r)
            if k is None:
                continue
            seen[r["CB"]][k].add(r["UB"])
    return {b: {k: len(v) for k, v in d.items()} for b, d in seen.items()}


def add(dst, ex):
    for bc, d in ex.items():
        t = dst.setdefault(bc, {})
        for k, v in d.items():
            t[k] = t.get(k, 0) + v


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    meta, mat = load_meta(), matrix_counts()
    samples = sorted(meta, key=lambda s: int("".join(c for c in s[3:5] if c.isdigit())))

    per_state = {}
    for state in ("base", "resc", "anti"):
        rows, pair_obs, pair_exp, region_obs = [], defaultdict(int), defaultdict(float), defaultdict(int)
        n_orf_dist = defaultdict(lambda: defaultdict(int))
        for sample in samples:
            sex, diet, fmap = meta[sample]
            bcs = barcodes(sample)
            n_cells = len(bcs)
            cells = {bcs[c - 1]: dict(v) for c, v in mat.get(sample, {}).items()}
            if state in ("resc", "anti"):
                add(cells, extra(RESCUE, sample, "_rescued_reads.tsv",
                                 lambda r: r["construct"].replace("_cassette", "")
                                 if r["in_called_cell"] == "1" and r["CB"] and r["UB"] else None))
            if state == "anti":
                add(cells, extra(ANTI, sample, "_antisense_reads.tsv",
                                 lambda r: None if r["contig"] == "AAV_shared"
                                 else r["contig"].replace("_cassette", "")))

            marg = {o: sum(1 for c in cells.values() if c.get(o, 0) >= 1) for o in ORFS}
            dist = defaultdict(int)
            for c in cells.values():
                hot = [o for o in ORFS if c.get(o, 0) >= 1]
                dist[len(hot)] += 1
                if len(hot) >= 2:
                    for a, b in combinations(sorted(hot), 2):
                        pair_obs[(a, b)] += 1
                        region_obs[tuple(sorted((fmap[a], fmap[b])))] += 1
                    if len(hot) == 3:
                        pair_obs[("ALL3", "")] += 1
            dist[0] += n_cells - len(cells)
            for k in (0, 1, 2, 3):
                n_orf_dist[sample][k] = dist[k]
            for a, b in combinations(ORFS, 2):
                key = tuple(sorted((a, b)))
                pair_exp[key] += n_cells * (marg[a] / n_cells) * (marg[b] / n_cells)
            rows.append((sample, sex, diet, n_cells, dist))
        per_state[state] = (rows, pair_obs, pair_exp, region_obs, n_orf_dist)

    rows, pair_obs, pair_exp, region_obs, n_orf_dist = per_state["anti"]

    print("Reporters per nucleus (final state: matrix + rescue + antisense)")
    print("=" * 76)
    print("%-11s %8s %10s %8s %8s %8s %9s"
          % ("sample", "nuclei", "0 (none)", "1 ORF", "2 ORFs", "3 ORFs", "%multi"))
    T = defaultdict(int)
    for sample, sex, diet, n, d in rows:
        multi = d[2] + d[3]
        pos = d[1] + multi
        print("%-11s %8d %10d %8d %8d %8d %8.2f%%"
              % (sample, n, d[0], d[1], d[2], d[3],
                 100 * multi / pos if pos else 0.0))
        for k in (0, 1, 2, 3):
            T[k] += d[k]
        T["n"] += n
    multi = T[2] + T[3]
    pos = T[1] + multi
    print("%-11s %8d %10d %8d %8d %8d %8.2f%%"
          % ("TOTAL", T["n"], T[0], T[1], T[2], T[3], 100 * multi / pos))
    print()
    print("So %d of %d reporter-positive nuclei carry more than one reporter (%.2f%%);"
          % (multi, pos, 100 * multi / pos))
    print("%d carry all three." % T[3])

    print()
    print("Is that more than coincidence? expected = n_cells * p_A * p_B, summed over samples")
    print("=" * 76)
    print("%-22s %9s %10s %8s" % ("reporter pair", "observed", "expected", "obs/exp"))
    for a, b in combinations(ORFS, 2):
        key = tuple(sorted((a, b)))
        o, e = pair_obs[key], pair_exp[key]
        print("%-22s %9d %10.1f %8s"
              % (a + "+" + b, o, e, "%.1f" % (o / e) if e else "n/a"))
    print()
    print("Region pairs (translated per mouse; the fluorophore map is randomised)")
    print("%-22s %9s" % ("region pair", "nuclei"))
    for k in sorted(region_obs, key=lambda x: -region_obs[x]):
        print("%-22s %9d" % ("+".join(k), region_obs[k]))

    print()
    print("Did the rescue tracks manufacture the overlaps?")
    print("=" * 76)
    print("%-22s %9s %9s %9s" % ("reporter pair", "base", "+resc", "+anti"))
    for a, b in combinations(ORFS, 2):
        key = tuple(sorted((a, b)))
        print("%-22s %9d %9d %9d"
              % (a + "+" + b, per_state["base"][1][key],
                 per_state["resc"][1][key], per_state["anti"][1][key]))

    with open(os.path.join(OUTDIR, "multiorf_per_sample.tsv"), "w") as fh:
        fh.write("sample\tsex\tdiet\ttotal_cells\tn_orf_0\tn_orf_1\tn_orf_2\tn_orf_3\n")
        for sample, sex, diet, n, d in rows:
            fh.write("%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\n"
                     % (sample, sex, diet, n, d[0], d[1], d[2], d[3]))
    with open(os.path.join(OUTDIR, "multiorf_pairs.tsv"), "w") as fh:
        fh.write("pair\tobserved\texpected_independent\tobs_over_exp\n")
        for a, b in combinations(ORFS, 2):
            key = tuple(sorted((a, b)))
            o, e = pair_obs[key], pair_exp[key]
            fh.write("%s+%s\t%d\t%.2f\t%s\n" % (a, b, o, e, "%.3f" % (o / e) if e else "NA"))
    print()
    print("written to", OUTDIR)


if __name__ == "__main__":
    main()
