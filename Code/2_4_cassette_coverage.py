#!/usr/bin/env python3
"""Step 2.4 stage B: map reads onto the reconstructed cassettes by exact 31-mer
match and report where they land.

Reads SAM records on stdin.  Usage:
    2_4_cassette_coverage.py <sample> <stream_label> <outdir>

The point is the split by `in_v1_reference`: a read whose every matched base
lies outside the v1 contigs is a read the current reference could not have
counted under any threshold, and one the 2_2 unmapped scan could not have seen
either - that scan only searched for k-mers of the v1 contigs themselves.
"""
import os, sys, collections

PROJ = "/scratch/user/uqschauq/Insular"
BLOCKS = os.path.join(PROJ, "Results/2_4_cassette_blocks")
K = 31

comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(comp)[::-1]


def load_cassettes():
    d, cur = {}, None
    for line in open(os.path.join(BLOCKS, "cassettes.fa")):
        if line.startswith(">"):
            cur = line[1:].split()[0]; d[cur] = []
        elif cur:
            d[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


def load_mask():
    """(cassette, pos) -> (in_v1, n_constructs).  pos is 1-based."""
    modelled, nconstr = {}, {}
    with open(os.path.join(BLOCKS, "cassette_basemask.tsv")) as fh:
        next(fh)
        for line in fh:
            c, p, m, n, _feat = line.rstrip("\n").split("\t")
            modelled[(c, int(p))] = int(m)
            nconstr[(c, int(p))] = int(n)
    return modelled, nconstr


def main():
    sample, stream, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)

    cassettes = load_cassettes()
    modelled, nconstr = load_mask()

    # 31-mer -> list of (cassette, start_pos_1based).  Both strands of the
    # cassette are indexed so a read from either orientation is found.
    index = collections.defaultdict(list)
    for name, s in cassettes.items():
        L = len(s)
        for i in range(L - K + 1):
            index[s[i:i + K]].append((name, i + 1))
        r = rc(s)
        for i in range(L - K + 1):
            index[r[i:i + K]].append((name, L - i - K + 1))
    index = dict(index)
    sys.stderr.write(f"{sample}/{stream}: indexed {len(index)} 31-mers "
                     f"over {len(cassettes)} cassettes\n")

    depth = collections.Counter()          # (cassette, pos) -> reads
    n = hit = 0
    cls = collections.Counter()            # read classification
    per_cassette = collections.Counter()

    for line in sys.stdin:
        if line[0] == "@":
            continue
        f = line.split("\t", 11)
        if len(f) < 10:
            continue
        n += 1
        s = f[9]
        pos = set()
        for i in range(len(s) - K + 1):
            m = index.get(s[i:i + K])
            if m:
                for (cname, p) in m:
                    pos.update((cname, q) for q in range(p, p + K))
        if not pos:
            continue
        hit += 1
        for key in pos:
            depth[key] += 1
        hitc = {c for (c, _) in pos}
        for c in hitc:
            per_cassette[c] += 1
        mod = [modelled.get(k, 0) for k in pos]
        if not any(mod):
            cls["entirely_outside_v1"] += 1
        elif not all(mod):
            cls["spans_v1_boundary"] += 1
        else:
            cls["entirely_inside_v1"] += 1
        # assignability: best (lowest) construct-count over the matched bases
        best = min(nconstr.get(k, 3) for k in pos)
        cls[f"assignable_to_{best}_construct" + ("" if best == 1 else "s")] += 1

    # ---- per-base coverage -------------------------------------------------
    cov = os.path.join(outdir, f"{sample}_{stream}_coverage.tsv")
    with open(cov, "w") as fh:
        fh.write("cassette\tpos\tdepth\tin_v1_reference\tn_constructs_k90\n")
        for name, s in cassettes.items():
            for p in range(1, len(s) + 1):
                d = depth.get((name, p), 0)
                if d:
                    fh.write(f"{name}\t{p}\t{d}\t{modelled.get((name,p),0)}\t"
                             f"{nconstr.get((name,p),0)}\n")

    # ---- summary -----------------------------------------------------------
    summ = os.path.join(outdir, f"{sample}_{stream}_summary.tsv")
    with open(summ, "w") as fh:
        fh.write("sample\tstream\tmetric\tvalue\n")
        def w(metric, value):
            fh.write(f"{sample}\t{stream}\t{metric}\t{value}\n")
        w("reads_scanned", n)
        w("reads_matching_cassette", hit)
        for k in ("entirely_inside_v1", "spans_v1_boundary", "entirely_outside_v1"):
            w(k, cls[k])
        for k in sorted(cls):
            if k.startswith("assignable_to_"):
                w(k, cls[k])
        for c in sorted(cassettes):
            w(f"reads_on_{c}", per_cassette[c])
    sys.stderr.write(f"{sample}/{stream}: scanned {n}, matched {hit}\n")


if __name__ == "__main__":
    main()
