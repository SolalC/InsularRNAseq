#!/usr/bin/env python3
"""
Step 2.8: join the 2.7 BLAST hits back to cell barcodes and UMIs.

The reads BLAST found among the unmapped pool still carry Cell Ranger's
corrected tags (verified on Ech1FSD: CB on 91.1% of unmapped reads, UB on
99.5%), so a rescued read can be turned back into a molecule in a called cell.

What is kept (the same rule 2.7's summary used, so the numbers line up):
  - confident HSPs only: >= MIN_LEN bp at >= MIN_PID identity
  - one HSP per read per construct (tdTomato is a tandem dimer and hits itself
    twice; keeping both would double count)
  - a read is assigned only when its best construct beats the next construct by
    MARGIN bits - EGFP and EBFP2 are 97% identical over the ORF, so ties stay
    ambiguous and are never split
  - the HSP must cover >= MIN_DIAG bases that are construct-diagnostic, i.e.
    n_constructs_k90 == 1 in the 2.4 base mask - the same rule the v2 reference
    masking used. A WPRE/polyA hit with no diagnostic base assigns nothing.
  - HSPs falling entirely inside EBFP2_cassette 715-873 are dropped: step 2.6
    showed that window is mouse Ig heavy-chain V sequence at 0 mismatches.

THE AMBIENT CONTROL. Rescued reads are reads STAR rejected, so they bypass the
quality bar every mouse gene passed. The question that decides whether they are
worth anything is whether they concentrate in called cells or are smeared over
all barcodes like ambient RNA. While streaming the BAM this also counts, over
EVERY unmapped read, what fraction of barcodes are called cells (f_bg). If the
rescued reads show the same fraction (f_resc ~ f_bg) the signal is ambient and
buys nothing; if f_resc >> f_bg it is cellular. Same reads, same stream, no
assumption about empty-droplet denominators.

Reads are relocated in the BAM by their position in the `samtools view '*'`
stream, which 2.7 recorded in the query id as r<N>. The read name is verified
against the query id, so a stream that came back in a different order is a hard
error rather than a silent mis-join.
"""

import argparse
import csv
import gzip
import os
import subprocess
import sys
from collections import defaultdict

PROJ = "/scratch/user/uqschauq/Insular"
HITDIR = os.path.join(PROJ, "Results", "2_7_unmapped_blast")
COUNTS = os.path.join(PROJ, "Data", "Counts_v2")
BASEMASK = os.path.join(PROJ, "Results", "2_4_cassette_blocks", "cassette_basemask.tsv")
OUTDIR = os.path.join(PROJ, "Results", "2_8_rescue_assign")

MIN_LEN, MIN_PID, MARGIN, MIN_DIAG = 30, 90.0, 5.0, 20
IGH = ("EBFP2_cassette", 715, 873)
CONSTRUCTS = ["EGFP_cassette", "tdTomato_cassette", "EBFP2_cassette"]
SHORT = {c: c.replace("_cassette", "") for c in CONSTRUCTS}


def load_diag():
    diag = defaultdict(set)
    with open(BASEMASK) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["n_constructs_k90"] == "1":
                diag[r["cassette"] + "_cassette"].add(int(r["pos"]))
    return diag


def pick_reads(sample, diag):
    """Read -> (construct, diagnostic bases). Applies every 2.7 filter."""
    path = os.path.join(HITDIR, sample + "_hits.tsv")
    best, span = defaultdict(dict), {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            q, sub = f[0], f[1]
            pid, ln = float(f[2]), int(f[3])
            a, b = int(f[8]), int(f[9])
            bits = float(f[11])
            if ln < MIN_LEN or pid < MIN_PID:
                continue
            if sub not in best[q] or bits > best[q][sub][0]:
                best[q][sub] = (bits,)
                span[(q, sub)] = (min(a, b), max(a, b))

    keep = {}
    for q, hits in best.items():
        ranked = sorted(hits.items(), key=lambda kv: -kv[1][0])
        top = ranked[0][0]
        second = ranked[1][1][0] if len(ranked) > 1 else 0.0
        if ranked[0][1][0] - second < MARGIN:
            continue
        lo, hi = span[(q, top)]
        if top == IGH[0] and lo >= IGH[1] and hi <= IGH[2]:
            continue
        ndiag = len(diag[top] & set(range(lo, hi + 1)))
        if ndiag < MIN_DIAG:
            continue
        nr = int(q[1:].split("_", 1)[0])
        name = q.split("_", 1)[1]
        keep[nr] = (name, top, ndiag)
    return keep


def called_cells(sample):
    p = os.path.join(COUNTS, sample, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")
    with gzip.open(p, "rt") as fh:
        bcs = [l.strip() for l in fh]
    # column index in the matrix is 1-based line number, which is how
    # 2_1_projection_validation.R identifies a cell.
    return set(bcs), {b: i + 1 for i, b in enumerate(bcs)}


def existing_positives(sample, col_of):
    """(barcode, construct_short) already >=1 UMI in the v2 matrix."""
    p = os.path.join(PROJ, "Results", "2_1_projection_validation_v2",
                     "per_cell_calls_ge1.tsv.gz")
    bc_of = {v: k for k, v in col_of.items()}
    pos = set()
    with gzip.open(p, "rt") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["sample"] != sample:
                continue
            bc = bc_of.get(int(r["col"]))
            if bc is None:
                continue
            for c in ("EGFP", "tdTomato", "EBFP2"):
                if int(r[c]) >= 1:
                    pos.add((bc, c))
    return pos


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    args = ap.parse_args()
    sample = args.sample

    os.makedirs(OUTDIR, exist_ok=True)
    diag = load_diag()
    keep = pick_reads(sample, diag)
    cells, col_of = called_cells(sample)
    prior = existing_positives(sample, col_of)
    print("%s: %d reads pass the 2.7 filters, %d called cells"
          % (sample, len(keep), len(cells)))

    bam = os.path.join(COUNTS, sample, "outs", "possorted_genome_bam.bam")
    proc = subprocess.Popen(["samtools", "view", "-@", "2", bam, "*"],
                            stdout=subprocess.PIPE, text=True, bufsize=1 << 20)

    rows = []
    n_total = n_withcb = n_bg_cell = 0
    for nr, line in enumerate(proc.stdout, start=1):
        n_total += 1
        cb = ub = None
        # Only parse tags on lines we need, plus the CB for the background rate.
        tab = line.rstrip("\n").split("\t")
        for t in tab[11:]:
            if t.startswith("CB:Z:"):
                cb = t[5:]
            elif t.startswith("UB:Z:"):
                ub = t[5:]
        if cb:
            n_withcb += 1
            if cb in cells:
                n_bg_cell += 1
        hit = keep.get(nr)
        if hit is None:
            continue
        name, construct, ndiag = hit
        if tab[0] != name:
            sys.exit("FATAL %s: record %d is %s, expected %s - the BAM stream "
                     "came back in a different order, join aborted"
                     % (sample, nr, tab[0], name))
        rows.append((nr, name, construct, ndiag, cb or "", ub or "",
                     "1" if (cb and cb in cells) else "0"))
    proc.stdout.close()
    if proc.wait() != 0:
        sys.exit("FATAL: samtools exited %d" % proc.returncode)

    if len(rows) != len(keep):
        sys.exit("FATAL %s: relocated %d of %d reads" % (sample, len(rows), len(keep)))

    with open(os.path.join(OUTDIR, sample + "_rescued_reads.tsv"), "w") as fh:
        fh.write("read_nr\tread_name\tconstruct\tdiag_bases\tCB\tUB\tin_called_cell\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    # Molecules = distinct (CB, UB) in a called cell. Newly positive = that cell
    # had no UMI of that construct in the v2 matrix.
    mols = defaultdict(set)
    for _, _, construct, _, cb, ub, incell in rows:
        if incell == "1" and cb and ub:
            mols[construct].add((cb, ub))
    newcells = defaultdict(set)
    oldcells = defaultdict(set)
    for construct, mset in mols.items():
        for cb, _ in mset:
            (newcells if (cb, SHORT[construct]) not in prior else oldcells)[construct].add(cb)

    f_bg = n_bg_cell / n_withcb if n_withcb else 0.0
    n_resc_cb = sum(1 for r in rows if r[4])
    n_resc_cell = sum(1 for r in rows if r[6] == "1")
    f_resc = n_resc_cell / n_resc_cb if n_resc_cb else 0.0

    with open(os.path.join(OUTDIR, sample + "_rescue_summary.tsv"), "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write("unmapped_reads\t%d\n" % n_total)
        fh.write("unmapped_reads_with_CB\t%d\n" % n_withcb)
        fh.write("background_frac_CB_in_called_cell\t%.6f\n" % f_bg)
        fh.write("rescued_reads\t%d\n" % len(rows))
        fh.write("rescued_reads_with_CB\t%d\n" % n_resc_cb)
        fh.write("rescued_reads_in_called_cell\t%d\n" % n_resc_cell)
        fh.write("rescued_frac_CB_in_called_cell\t%.6f\n" % f_resc)
        fh.write("enrichment_over_background\t%.3f\n"
                 % (f_resc / f_bg if f_bg else float("nan")))
        for c in CONSTRUCTS:
            fh.write("molecules_%s\t%d\n" % (SHORT[c], len(mols.get(c, ()))))
            fh.write("cells_new_%s\t%d\n" % (SHORT[c], len(newcells.get(c, ()))))
            fh.write("cells_already_positive_%s\t%d\n" % (SHORT[c], len(oldcells.get(c, ()))))

    print("%s: %d rescued reads, %d in called cells (f_resc=%.4f vs f_bg=%.4f, "
          "enrichment %.2fx)" % (sample, len(rows), n_resc_cell, f_resc, f_bg,
                                 f_resc / f_bg if f_bg else float("nan")))
    for c in CONSTRUCTS:
        print("    %-9s molecules %4d  new cells %4d  already-positive cells %4d"
              % (SHORT[c], len(mols.get(c, ())), len(newcells.get(c, ())),
                 len(oldcells.get(c, ()))))


if __name__ == "__main__":
    main()
