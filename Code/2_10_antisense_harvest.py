#!/usr/bin/env python3
"""
Step 2.10: harvest the ANTISENSE reads already aligned to the transgene contigs.

The v2 GTF declares one `+`-strand exon per transgene contig, so a read aligning
antisense to a cassette gets Cell Ranger's `AN` tag instead of `GX` and never
reaches the count matrix. Those reads are already mapped, already carry a
corrected CB and UB, and in Ech3MSD amount to 694 of 10,577 transgene-aligned
reads (6.6%) - comparable to everything the 2.7 BLAST rescue recovered for the
same animal.

Antisense in 10x data is part genuine (ITR-driven antisense transcription is well
described for AAV) and part template-switching artefact. Either way it is
proportional to real sense expression, which is what a presence/absence call
needs. Whether to believe it is settled by the controls below, not by assertion.

Reads are pulled by indexed region query on the four small transgene contigs, so
this touches a few thousand records rather than streaming the 29 GB BAM.

Two controls, both computed in the same pass:

  1. MATCHED AMBIENT CONTROL. The fraction of ANTISENSE reads whose barcode is a
     called cell is compared with the same fraction for SENSE reads on the very
     same contigs. Sense reads are the signal we already trust and count, so if
     antisense lands in cells at the same rate it is behaving like real transgene
     signal rather than ambient. This is a tighter comparison than the
     unmapped-pool background used in 2.8, because both arms come from the same
     molecules on the same contigs.
  2. FAILED-REPORTER CONTROL (applied at aggregation, step 2.11). tdTomato failed
     biologically in Ech2/6/10/11/12. Antisense tdTomato molecules must stay near
     zero there, as the BLAST rescue did (4 vs 573).

AAV_shared antisense counts toward "transgene-positive" but assigns no
projection - the shared element carries no construct information by design.
"""

import argparse
import csv
import gzip
import os
import subprocess
import sys
from collections import defaultdict

PROJ = "/scratch/user/uqschauq/Insular"
COUNTS = os.path.join(PROJ, "Data", "Counts_v2")
OUTDIR = os.path.join(PROJ, "Results", "2_10_antisense_harvest")

CONTIGS = ["EGFP_cassette", "tdTomato_cassette", "EBFP2_cassette", "AAV_shared"]
ORF_CONTIGS = CONTIGS[:3]


def called_cells(sample):
    p = os.path.join(COUNTS, sample, "outs", "filtered_feature_bc_matrix",
                     "barcodes.tsv.gz")
    with gzip.open(p, "rt") as fh:
        return set(l.strip() for l in fh)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    args = ap.parse_args()
    sample = args.sample

    os.makedirs(OUTDIR, exist_ok=True)
    cells = called_cells(sample)
    bam = os.path.join(COUNTS, sample, "outs", "possorted_genome_bam.bam")
    if not os.path.exists(bam):
        sys.exit("FATAL: no BAM at " + bam)

    proc = subprocess.Popen(["samtools", "view", bam] + CONTIGS,
                            stdout=subprocess.PIPE, text=True, bufsize=1 << 20)

    rows = []
    n = defaultdict(int)
    for line in proc.stdout:
        f = line.rstrip("\n").split("\t")
        contig = f[2]
        if contig not in CONTIGS:
            continue
        cb = ub = None
        has_gx = has_an = False
        for t in f[11:]:
            if t.startswith("GX:Z:"):
                has_gx = True
            elif t.startswith("AN:Z:"):
                has_an = True
            elif t.startswith("CB:Z:"):
                cb = t[5:]
            elif t.startswith("UB:Z:"):
                ub = t[5:]
        sense = has_gx
        anti = has_an and not has_gx
        if not (sense or anti):
            n["neither"] += 1
            continue
        arm = "sense" if sense else "anti"
        n[arm] += 1
        n[arm + "_contig_" + contig] += 1
        if cb:
            n[arm + "_cb"] += 1
            if cb in cells:
                n[arm + "_cb_incell"] += 1
        if anti and cb and ub and cb in cells:
            rows.append((contig, cb, ub))
    proc.stdout.close()
    if proc.wait() != 0:
        sys.exit("FATAL: samtools exited %d" % proc.returncode)

    with open(os.path.join(OUTDIR, sample + "_antisense_reads.tsv"), "w") as fh:
        fh.write("contig\tCB\tUB\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")

    mols = defaultdict(set)
    for contig, cb, ub in rows:
        mols[contig].add((cb, ub))

    f_sense = n["sense_cb_incell"] / n["sense_cb"] if n["sense_cb"] else 0.0
    f_anti = n["anti_cb_incell"] / n["anti_cb"] if n["anti_cb"] else 0.0

    with open(os.path.join(OUTDIR, sample + "_antisense_summary.tsv"), "w") as fh:
        fh.write("metric\tvalue\n")
        fh.write("transgene_aligned_reads\t%d\n" % (n["sense"] + n["anti"] + n["neither"]))
        fh.write("sense_reads\t%d\n" % n["sense"])
        fh.write("antisense_reads\t%d\n" % n["anti"])
        fh.write("untagged_reads\t%d\n" % n["neither"])
        fh.write("antisense_frac\t%.6f\n"
                 % (n["anti"] / (n["sense"] + n["anti"]) if (n["sense"] + n["anti"]) else 0.0))
        fh.write("sense_frac_CB_in_called_cell\t%.6f\n" % f_sense)
        fh.write("antisense_frac_CB_in_called_cell\t%.6f\n" % f_anti)
        fh.write("anti_over_sense_incell_ratio\t%.4f\n"
                 % (f_anti / f_sense if f_sense else float("nan")))
        for c in CONTIGS:
            fh.write("antisense_reads_%s\t%d\n" % (c, n["anti_contig_" + c]))
            fh.write("antisense_molecules_%s\t%d\n" % (c, len(mols.get(c, ()))))

    print("%s: %d transgene-aligned reads, %d sense, %d antisense (%.1f%%), "
          "%d untagged" % (sample, n["sense"] + n["anti"] + n["neither"],
                           n["sense"], n["anti"],
                           100 * n["anti"] / (n["sense"] + n["anti"]) if (n["sense"] + n["anti"]) else 0.0,
                           n["neither"]))
    print("    in called cells: sense %.3f  antisense %.3f  ratio %.2f"
          % (f_sense, f_anti, f_anti / f_sense if f_sense else float("nan")))
    for c in CONTIGS:
        print("    %-18s antisense reads %4d  molecules in cells %4d"
              % (c, n["anti_contig_" + c], len(mols.get(c, ()))))


if __name__ == "__main__":
    main()
