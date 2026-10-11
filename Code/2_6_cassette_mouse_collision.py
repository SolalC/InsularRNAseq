#!/usr/bin/env python3
"""
Step 2.6: can cassette sequence collide with the mouse genome?

The v2 GTF annotates each transgene contig as a single `exon 1..<length>`, so a
UMI anywhere on a cassette - promoter, intron, ITR, WPRE - is counted as a
reporter UMI. At the 1-UMI presence/absence threshold that makes any
mouse<->cassette sequence homology a false-positive source, not merely a
read-loss source. The v2 masking only enforced uniqueness AMONG THE THREE
CONSTRUCTS; the mouse genome was never a competitor in that calculation.

The cassettes are not purely synthetic. From the plasmid maps: hCMV enhancer,
chicken beta-actin promoter + intron, rabbit beta-globin intron 1 and exon 3
(the CAG promoter of p24/p131), the human Synapsin 1 promoter (p140/EBFP2) and
Promega's chimeric intron (human beta-globin/IgG). Those are vertebrate and
could in principle align to mouse; the AAV ITR, WPRE and SV40 polyA are viral
and cannot.

Method - `tile`:
  Every 90 bp window (= the R2 read length) of each construct, step 1 bp, is
  written as a pseudo-read. The three reporter cassettes are tiled from the
  UNMASKED sequences in Results/2_4_cassette_blocks/cassettes.fa, because that
  is what a real read off the construct looks like; the N-masking exists only
  in the reference. AAV_shared is tiled from its reference contig with any
  N-containing window dropped (its spacers are artificial).

  Each tile is annotated with the fraction of its 90 positions that survive
  unmasked in the v2 reference contig ("kept"), and with the plasmid features
  it overlaps.

Method - `report`:
  Tiles are aligned by the caller against the MOUSE-ONLY index
  (refdata-gex-GRCm39-2024-A/star) under deliberately permissive filters, a
  superset of whatever Cell Ranger uses internally, so a clean result closes the
  question regardless of Cell Ranger's exact parameters. This step parses the
  SAM and reports, per tile, the best mouse hit and its mismatch count.

Interpreting the result - why one alignment run answers both directions:

  A real construct read matches its cassette with 0 mismatches; a real mouse
  read matches its own locus with 0 mismatches. So a read can only be diverted -
  in either direction - where cassette and mouse are essentially identical over
  the 90 bp. A tile whose best mouse hit carries several mismatches cannot steal
  a construct read (it loses to the perfect cassette hit) and cannot manufacture
  a false reporter UMI from a mouse read (it loses to the perfect mouse hit).
  Tiles with <= MM_RISK mismatches are therefore flagged; everything else is
  noise. Windows that are masked in v2 are reported separately: a collision
  there costs reads but cannot create a false projection call, because those
  bases are unmappable on the cassette anyway.
"""

import argparse
import os
import sys
from collections import defaultdict

PROJ = "/scratch/user/uqschauq/Insular"
CASSETTES = os.path.join(PROJ, "Results", "2_4_cassette_blocks", "cassettes.fa")
FEATURES = os.path.join(PROJ, "Results", "2_4_cassette_blocks", "features.tsv")
V2FASTA = os.path.join(PROJ, "Data", "Reference", "GRCm39_cassette_v2", "fasta", "genome.fa")

READLEN = 90          # R2 length; a tile is exactly one possible read
STEP = 1
MM_RISK = 5           # <= this many mismatches over 90 bp = a real risk
CONTIGS = ["EGFP_cassette", "tdTomato_cassette", "EBFP2_cassette", "AAV_shared"]


def read_fasta(path, keep=None):
    seqs, name, buf = {}, None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None and (keep is None or name in keep):
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name is not None and (keep is None or name in keep):
        seqs[name] = "".join(buf)
    return seqs


def load_features():
    """position -> set of plasmid feature names, per cassette short name."""
    feat = defaultdict(lambda: defaultdict(set))
    with open(FEATURES) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            cas = f[idx["cassette"]]
            name = f[idx["feature"]]
            for p in range(int(f[idx["start"]]), int(f[idx["end"]]) + 1):
                feat[cas][p].add(name)
    return feat


def cmd_tile(outdir):
    os.makedirs(outdir, exist_ok=True)

    # Unmasked construct sequence for the three reporters.
    cassettes = read_fasta(CASSETTES)
    # v2 reference contigs, to learn which bases survived masking.
    v2 = read_fasta(V2FASTA, keep=set(CONTIGS))
    missing = [c for c in CONTIGS if c not in v2]
    if missing:
        sys.exit("v2 reference is missing contig(s): " + ", ".join(missing))

    feat = load_features()

    fq_path = os.path.join(outdir, "tiles.fq")
    an_path = os.path.join(outdir, "tiles_annot.tsv")
    n_written = n_skipped = 0

    with open(fq_path, "w") as fq, open(an_path, "w") as an:
        an.write("tile\tcontig\tstart\tend\tkept_bases\tfeatures\n")
        for contig in CONTIGS:
            short = contig.replace("_cassette", "")
            # Reporters: tile the real construct. AAV_shared has no unmasked
            # counterpart (it is assembled from p24 blocks), so tile the contig.
            # cassettes.fa is keyed by the full contig name; features.tsv by the short one.
            src = cassettes.get(contig) if contig.endswith("_cassette") else v2[contig]
            if src is None:
                sys.exit("no unmasked sequence for " + contig)
            if contig.endswith("_cassette") and len(src) != len(v2[contig]):
                sys.exit("%s: unmasked %d bp vs v2 contig %d bp - not comparable"
                         % (contig, len(src), len(v2[contig])))

            ref = v2[contig]
            for start in range(0, len(src) - READLEN + 1, STEP):
                tile = src[start:start + READLEN]
                if "N" in tile:          # artificial spacer in AAV_shared
                    n_skipped += 1
                    continue
                kept = sum(1 for b in ref[start:start + READLEN] if b != "N")
                names = set()
                for p in range(start + 1, start + READLEN + 1):
                    names |= feat[short].get(p, set())
                tid = "%s:%d" % (contig, start + 1)
                fq.write("@%s\n%s\n+\n%s\n" % (tid, tile, "I" * READLEN))
                an.write("%s\t%s\t%d\t%d\t%d\t%s\n" % (
                    tid, contig, start + 1, start + READLEN, kept,
                    ",".join(sorted(names)) if names else "-"))
                n_written += 1

    print("tiles written : %d  (%s)" % (n_written, fq_path))
    print("tiles skipped : %d  (contained N)" % n_skipped)


def cmd_report(outdir, sam):
    annot = {}
    with open(os.path.join(outdir, "tiles_annot.tsv")) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            annot[f[0]] = dict(contig=f[1], start=int(f[2]), end=int(f[3]),
                               kept=int(f[4]), features=f[5])

    # Best (fewest-mismatch) mouse hit per tile.
    best = {}
    with open(sam) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            f = line.rstrip("\n").split("\t")
            qname, flag, rname, pos = f[0], int(f[1]), f[2], f[3]
            if flag & 0x4 or rname == "*":
                continue
            nm = None
            for tag in f[11:]:
                if tag.startswith("nM:i:"):      # STAR's per-pair mismatch tag
                    nm = int(tag[5:])
                elif tag.startswith("NM:i:") and nm is None:
                    nm = int(tag[5:])
            if nm is None:
                nm = 99
            prev = best.get(qname)
            if prev is None or nm < prev[0]:
                best[qname] = (nm, rname, pos)

    total = len(annot)
    aligned = len(best)
    print("=" * 72)
    print("Cassette 90-mers aligned against the MOUSE-ONLY index")
    print("=" * 72)
    print("tiles tested      : %d" % total)
    print("tiles with any mouse alignment : %d (%.4f%%)"
          % (aligned, 100.0 * aligned / total if total else 0.0))
    print()

    if not aligned:
        print("No cassette window aligns anywhere in GRCm39 under permissive")
        print("filters. Neither direction of read diversion is possible:")
        print("  - no construct read can be stolen by a mouse locus")
        print("  - no mouse read can be converted into a reporter UMI")
        print("The transgene features and the mouse gene counts are independent.")
        return

    # Per-contig breakdown, split by whether the window is mappable in v2.
    print("--- mouse-aligning tiles per contig ---")
    print("%-18s %8s %8s %8s %8s" % ("contig", "tiles", "aligned", "mm<=%d" % MM_RISK, "risk_kept"))
    per = defaultdict(lambda: [0, 0, 0, 0])
    for tid, a in annot.items():
        per[a["contig"]][0] += 1
        if tid in best:
            per[a["contig"]][1] += 1
            if best[tid][0] <= MM_RISK:
                per[a["contig"]][2] += 1
                # "kept" = the window is still mappable on the cassette in v2,
                # so a collision there can create a false reporter UMI.
                if a["kept"] >= READLEN * 0.9:
                    per[a["contig"]][3] += 1
    for contig in CONTIGS:
        t, al, mm, risk = per[contig]
        print("%-18s %8d %8d %8d %8d" % (contig, t, al, mm, risk))
    print()

    risky = [(best[t][0], t, annot[t], best[t]) for t in best
             if best[t][0] <= MM_RISK]
    risky.sort()
    print("--- tiles within %d mismatches of a mouse locus (%d) ---"
          % (MM_RISK, len(risky)))
    if not risky:
        print("none - every mouse hit is far enough away that the true source")
        print("locus always wins the alignment in both directions.")
    else:
        print("%-28s %6s %6s %-10s %-12s %s"
              % ("tile", "mm", "kept", "mouse", "pos", "features"))
        for nm, tid, a, (mm, rname, pos) in risky[:200]:
            print("%-28s %6d %6d %-10s %-12s %s"
                  % (tid, mm, a["kept"], rname, pos, a["features"]))
        if len(risky) > 200:
            print("... %d more, see the full table" % (len(risky) - 200))

    # Full table for anything that aligned at all.
    out = os.path.join(outdir, "mouse_hits.tsv")
    with open(out, "w") as fh:
        fh.write("tile\tcontig\tstart\tend\tkept_bases\tmismatches\t"
                 "mouse_chrom\tmouse_pos\tfeatures\n")
        for tid in sorted(best, key=lambda k: (annot[k]["contig"], annot[k]["start"])):
            a = annot[tid]
            nm, rname, pos = best[tid]
            fh.write("%s\t%s\t%d\t%d\t%d\t%d\t%s\t%s\t%s\n"
                     % (tid, a["contig"], a["start"], a["end"], a["kept"],
                        nm, rname, pos, a["features"]))
    print("\nfull hit table -> %s" % out)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    t = sub.add_parser("tile");   t.add_argument("--outdir", required=True)
    r = sub.add_parser("report"); r.add_argument("--outdir", required=True)
    r.add_argument("--sam", required=True)
    args = ap.parse_args()
    if args.cmd == "tile":
        cmd_tile(args.outdir)
    else:
        cmd_report(args.outdir, args.sam)
