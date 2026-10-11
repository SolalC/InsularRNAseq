#!/usr/bin/env python3
"""
Step 2.7 (summary): fold the per-sample BLAST hit tables into a verdict.

Reads Results/2_7_unmapped_blast/<sample>_hits.tsv (blastn outfmt 6 against the
three full ITR-to-ITR cassettes) and answers: does BLAST find transgene sequence
among the unmapped reads that the exact-31-mer scans of 2_2 and 2_4 missed?

Three things the smoke test showed have to be handled, or the counts inflate:

  1. tdTomato is a TANDEM dimer. An ORF read hits its own cassette twice (the
     second unit at ~86% identity). Collapse to one HSP per read per construct.
  2. EGFP and EBFP2 are 97% identical over the ORF, so an EBFP2 read also hits
     EGFP_cassette at ~95%. A read is only assigned to a construct when its best
     hit beats the best hit to any OTHER construct by MARGIN bits; otherwise it
     is counted as ambiguous, never split between the two.
  3. WPRE and the polyA are near-identical across constructs by design, so most
     hits are expected to be ambiguous. That is not a failure - it is the same
     reason AAV_shared exists in the reference.

Tiers are reported rather than one cutoff, because the whole point of using
BLAST is to see what sits below the exact-match threshold.
"""

import csv
import glob
import os
import sys
from collections import defaultdict

PROJ = "/scratch/user/uqschauq/Insular"
HITDIR = os.path.join(PROJ, "Results", "2_7_unmapped_blast")
FEATURES = os.path.join(PROJ, "Results", "2_4_cassette_blocks", "features.tsv")

MARGIN = 5.0          # bits a construct must win by to own a read
CONSTRUCTS = ["EGFP_cassette", "tdTomato_cassette", "EBFP2_cassette"]

# The window of EBFP2_cassette that is mouse Ig heavy-chain V sequence at 0
# mismatches (step 2.6). A hit here is not evidence of transgene.
IGH_WINDOW = ("EBFP2_cassette", 715, 873)

TIERS = [
    ("permissive  (e<=1e-5)",          lambda L, p: True),
    ("confident   (>=30bp, >=90%)",    lambda L, p: L >= 30 and p >= 90.0),
    ("strict      (>=50bp, >=95%)",    lambda L, p: L >= 50 and p >= 95.0),
    ("near-exact  (>=70bp, >=99%)",    lambda L, p: L >= 70 and p >= 99.0),
]


def load_features():
    feat = defaultdict(lambda: defaultdict(set))
    with open(FEATURES) as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        i = {c: n for n, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cas = f[i["cassette"]] + "_cassette"
            for p in range(int(f[i["start"]]), int(f[i["end"]]) + 1):
                feat[cas][p].add(f[i["feature"]])
    return feat


def summarise(path, feat):
    """Return per-tier stats for one sample."""
    # best bitscore per (read, construct)
    best = defaultdict(dict)
    span = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, s = f[0], f[1]
            pid, ln = float(f[2]), int(f[3])
            sstart, send = int(f[8]), int(f[9])
            bits = float(f[11])
            prev = best[q].get(s)
            if prev is None or bits > prev[0]:
                best[q][s] = (bits, pid, ln)
                span[(q, s)] = (min(sstart, send), max(sstart, send))

    out = {}
    for label, keep in TIERS:
        assigned = defaultdict(int)
        ambiguous = 0
        igh = 0
        featcount = defaultdict(int)
        for q, hits in best.items():
            ok = {s: v for s, v in hits.items() if keep(v[2], v[1])}
            if not ok:
                continue
            ranked = sorted(ok.items(), key=lambda kv: -kv[1][0])
            top, topv = ranked[0]
            second = ranked[1][1][0] if len(ranked) > 1 else 0.0
            if topv[0] - second < MARGIN:
                ambiguous += 1
                continue
            lo, hi = span[(q, top)]
            if top == IGH_WINDOW[0] and lo >= IGH_WINDOW[1] and hi <= IGH_WINDOW[2]:
                igh += 1
                continue
            assigned[top] += 1
            names = set()
            for p in range(lo, hi + 1):
                names |= feat[top].get(p, set())
            featcount[",".join(sorted(names)) if names else "unannotated"] += 1
        out[label] = dict(assigned=dict(assigned), ambiguous=ambiguous,
                          igh=igh, features=dict(featcount))
    return out, len(best)


def main():
    feat = load_features()
    paths = sorted(glob.glob(os.path.join(HITDIR, "*_hits.tsv")))
    if not paths:
        sys.exit("no hit tables in " + HITDIR)

    per_sample = {}
    for p in paths:
        s = os.path.basename(p).replace("_hits.tsv", "")
        per_sample[s], nreads = summarise(p, feat)
        print("%-11s %8d reads with >=1 HSP" % (s, nreads))
    print()

    for label, _ in TIERS:
        print("=" * 74)
        print("TIER %s   (margin %.0f bits to own a read)" % (label, MARGIN))
        print("=" * 74)
        print("%-11s %9s %9s %9s %11s %9s"
              % ("sample", "EGFP", "tdTomato", "EBFP2", "ambiguous", "IgH_win"))
        tot = defaultdict(int)
        for s in sorted(per_sample):
            d = per_sample[s][label]
            a = d["assigned"]
            print("%-11s %9d %9d %9d %11d %9d"
                  % (s, a.get("EGFP_cassette", 0), a.get("tdTomato_cassette", 0),
                     a.get("EBFP2_cassette", 0), d["ambiguous"], d["igh"]))
            for k in CONSTRUCTS:
                tot[k] += a.get(k, 0)
            tot["ambiguous"] += d["ambiguous"]
            tot["igh"] += d["igh"]
        print("%-11s %9d %9d %9d %11d %9d"
              % ("TOTAL", tot["EGFP_cassette"], tot["tdTomato_cassette"],
                 tot["EBFP2_cassette"], tot["ambiguous"], tot["igh"]))
        print()

    # Where on the cassette the confident hits land - promoter/intron hits are
    # the ones that would not have been diagnostic anyway.
    print("=" * 74)
    print("Cassette region of confident (>=30bp, >=90%) assigned hits, pooled")
    print("=" * 74)
    pooled = defaultdict(int)
    for s in per_sample:
        for k, v in per_sample[s]["confident   (>=30bp, >=90%)"]["features"].items():
            pooled[k] += v
    for k, v in sorted(pooled.items(), key=lambda kv: -kv[1]):
        print("%8d  %s" % (v, k))


if __name__ == "__main__":
    main()
