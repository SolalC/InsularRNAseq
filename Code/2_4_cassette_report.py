#!/usr/bin/env python3
"""Step 2.4 report: aggregate stage B coverage over the reconstructed cassettes.

Recomputes v1 mappability properly. Stage A's `in_v1_reference` column marks a
base only when a whole v1 contig matches exactly, which understates reality --
EBFP2's WPRE shares a 589 bp prefix with the p24-derived shared contig, so reads
there do map today. Here a base counts as v1-mappable if the read-length window
starting on it exists in the v1 contigs at all.

  k=90 : strict, read maps end to end
  k=60 : permissive, ~the shortest match that survives Cell Ranger's
         --outFilterMatchNminOverLread 0.66 on a 90 bp read via soft-clipping
"""
import os, collections, glob

PROJ = "/scratch/user/uqschauq/Insular"
BLOCKS = os.path.join(PROJ, "Results/2_4_cassette_blocks")
COV = os.path.join(PROJ, "Results/2_4_cassette_coverage")
OUT = os.path.join(PROJ, "Results/2_4_cassette_report")
V1FA = os.path.join(PROJ, "Data/Reference/transgenes/transgenes.fa")

comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(comp)[::-1]

def load_fasta(p):
    d, cur = {}, None
    for line in open(p):
        if line.startswith(">"):
            cur = line[1:].split()[0]; d[cur] = []
        elif cur:
            d[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in d.items()}

def kmers(seqs, k):
    out = set()
    for s in seqs:
        for st in (s, rc(s)):
            for i in range(len(st) - k + 1):
                out.add(st[i:i + k])
    return out

def mappable_mask(cassettes, v1seqs, k):
    ks = kmers(v1seqs.values(), k)
    m = {}
    for c, s in cassettes.items():
        lab = [(s[i:i + k] in ks) for i in range(len(s) - k + 1)]
        lab += [lab[-1]] * (k - 1) if lab else []
        m[c] = lab
    return m

def main():
    os.makedirs(OUT, exist_ok=True)
    cassettes = load_fasta(os.path.join(BLOCKS, "cassettes.fa"))
    cassettes = {k.replace("_cassette", ""): v for k, v in cassettes.items()}
    v1 = load_fasta(V1FA)

    map90 = mappable_mask(cassettes, v1, 90)
    map60 = mappable_mask(cassettes, v1, 60)

    feat, nconstr = {}, {}
    with open(os.path.join(BLOCKS, "cassette_basemask.tsv")) as fh:
        next(fh)
        for line in fh:
            c, p, _m, n, f = line.rstrip("\n").split("\t")
            feat[(c, int(p))] = f
            nconstr[(c, int(p))] = int(n)

    depth = collections.Counter()
    samples = set()
    for path in sorted(glob.glob(os.path.join(COV, "*_coverage.tsv"))):
        base = os.path.basename(path)
        stream = "unmapped" if base.endswith("_unmapped_coverage.tsv") else "aligned"
        samples.add(base.split("_" + stream)[0])
        with open(path) as fh:
            next(fh)
            for line in fh:
                c, p, d, _v, _n = line.rstrip("\n").split("\t")
                depth[(stream, c.replace("_cassette", ""), int(p))] += int(d)

    lines = []
    def w(*a):
        s = " ".join(str(x) for x in a)
        print(s); lines.append(s)

    w(f"samples aggregated: {len(samples)}")
    w("")
    w("=== 1. is the cassette sequence the reads land on present in the v1 reference? ===")
    w("read-bases = sum over reads of the cassette bases each read matched")
    w(f"{'stream':<9} {'cassette':<9} {'read_bases':>12} {'v1_mappable':>12} "
      f"{'NOT_mappable':>13} {'pct_lost':>9}")
    for stream in ("aligned", "unmapped"):
        for c in sorted(cassettes):
            tot = ok = 0
            for p in range(1, len(cassettes[c]) + 1):
                d = depth.get((stream, c, p), 0)
                if not d: continue
                tot += d
                if map90[c][p - 1]: ok += d
            if tot:
                w(f"{stream:<9} {c:<9} {tot:>12} {ok:>12} {tot-ok:>13} "
                  f"{100.0*(tot-ok)/tot:>8.1f}%")

    w("")
    w("=== 2. where on the cassette do the reads actually land? ===")
    w("(both streams pooled; this is the capture-geometry answer)")
    order = ["AAV-2_ITR", "Enhancer", "enhancer", "promoter", "intron", "exon",
             "Kozak", "ORF", "WPRE", "polyA_signal", "-"]
    for c in sorted(cassettes):
        agg = collections.Counter()
        for p in range(1, len(cassettes[c]) + 1):
            d = depth.get(("aligned", c, p), 0) + depth.get(("unmapped", c, p), 0)
            if d: agg[feat.get((c, p), "-")] += d
        tot = sum(agg.values())
        if not tot: continue
        w(f"  -- {c}_cassette   total read-bases {tot}")
        for f in order:
            if agg[f]:
                w(f"     {f:<14} {agg[f]:>12}  {100.0*agg[f]/tot:>6.1f}%")

    w("")
    w("=== 3. how much of the signal is construct-assignable? ===")
    w("n_constructs = how many cassettes share the 90 bp window at that base")
    w(f"{'cassette':<9} {'unique(1)':>13} {'2-way':>12} {'3-way':>12} {'pct_unique':>11}")
    for c in sorted(cassettes):
        agg = collections.Counter()
        for p in range(1, len(cassettes[c]) + 1):
            d = depth.get(("aligned", c, p), 0) + depth.get(("unmapped", c, p), 0)
            if d: agg[nconstr.get((c, p), 0)] += d
        tot = sum(agg.values())
        if tot:
            w(f"{c:<9} {agg[1]:>13} {agg[2]:>12} {agg[3]:>12} "
              f"{100.0*agg[1]/tot:>10.1f}%")

    w("")
    w("=== 4. the recoverable pool ===")
    w("read-bases on sequence that is BOTH unique to one construct AND absent")
    w("from the v1 reference: projection-informative signal currently discarded.")
    w(f"{'cassette':<9} {'strict k=90':>13} {'permissive k=60':>17}  top regions")
    for c in sorted(cassettes):
        r90 = r60 = 0
        byfeat = collections.Counter()
        for p in range(1, len(cassettes[c]) + 1):
            d = depth.get(("aligned", c, p), 0) + depth.get(("unmapped", c, p), 0)
            if not d or nconstr.get((c, p), 0) != 1: continue
            if not map90[c][p - 1]:
                r90 += d; byfeat[feat.get((c, p), "-")] += d
            if not map60[c][p - 1]:
                r60 += d
        top = ", ".join(f"{f}:{n}" for f, n in byfeat.most_common(3))
        w(f"{c:<9} {r90:>13} {r60:>17}  {top}")

    with open(os.path.join(OUT, "report.txt"), "w") as fh:
        fh.write("\n".join(lines) + "\n")

    with open(os.path.join(OUT, "per_base.tsv"), "w") as fh:
        fh.write("cassette\tpos\tfeature\tn_constructs\tv1_mappable_k90\t"
                 "v1_mappable_k60\tdepth_aligned\tdepth_unmapped\n")
        for c in sorted(cassettes):
            for p in range(1, len(cassettes[c]) + 1):
                fh.write(f"{c}\t{p}\t{feat.get((c,p),'-')}\t{nconstr.get((c,p),0)}\t"
                         f"{int(map90[c][p-1])}\t{int(map60[c][p-1])}\t"
                         f"{depth.get(('aligned',c,p),0)}\t{depth.get(('unmapped',c,p),0)}\n")
    print(f"\nwritten to {OUT}")

if __name__ == "__main__":
    main()
