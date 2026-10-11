#!/usr/bin/env python3
"""Step 2.4 stage A: reconstruct the packaged AAV cassette for each construct
from the plasmid maps, then partition the three cassettes into blocks by
sequence sharing.

Inputs : Data/Plasmid/{p24,p131,p140}.gb   (the actual constructs used)
Outputs: Results/2_4_cassette_blocks/

Why: the v1 reference models four hand-cut fragments (three ORFs + one shared
3' element taken from p24). Anything in the cassette outside those fragments is
unmappable and invisible to the 2_2 unmapped scan, which only searched for
k-mers of the fragments themselves. This measures what was left out, and which
cassette regions would multimap if full-length cassettes were used instead.
"""
import os, re, sys, collections

PROJ = "/scratch/user/uqschauq/Insular"
PLASMID = os.path.join(PROJ, "Data/Plasmid")
V1FA = os.path.join(PROJ, "Data/Reference/transgenes/transgenes.fa")
OUT = os.path.join(PROJ, "Results/2_4_cassette_blocks")

# plasmid -> reporter carried by its ORF
CONSTRUCTS = [("EGFP", "p24.gb"), ("tdTomato", "p131.gb"), ("EBFP2", "p140.gb")]

K_INDEX = 31   # k-mer size for the read-scanning index handed to stage B
K_READ  = 90   # R2 length: a window this size shared between cassettes multimaps

comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def rc(s): return s.translate(comp)[::-1]


def parse_genbank(path):
    """Minimal GenBank reader: returns (length, [(key, start, end, strand)], sequence).
    Coordinates are returned 1-based inclusive, as written in the file."""
    seq_parts, feats, in_feat, in_origin = [], [], False, False
    length = None
    for line in open(path):
        if line.startswith("LOCUS"):
            m = re.search(r"\s(\d+)\s+bp", line)
            if m:
                length = int(m.group(1))
        elif line.startswith("FEATURES"):
            in_feat = True
        elif line.startswith("ORIGIN"):
            in_feat, in_origin = False, True
        elif in_origin:
            if line.startswith("//"):
                break
            seq_parts.append(re.sub(r"[^A-Za-z]", "", line))
        elif in_feat:
            # feature lines start with exactly 5 spaces then a key; qualifiers are deeper
            m = re.match(r"^ {5}(\S+)\s+(\S.*)$", line.rstrip("\n"))
            if not m:
                continue
            key, loc = m.group(1), m.group(2).strip()
            strand = "-" if loc.startswith("complement(") else "+"
            m2 = re.search(r"(\d+)\.\.(\d+)", loc)
            if not m2:
                continue
            feats.append((key, int(m2.group(1)), int(m2.group(2)), strand))
    seq = "".join(seq_parts).upper()
    if length is not None and len(seq) != length:
        sys.stderr.write(f"WARN {os.path.basename(path)}: LOCUS says {length} bp, "
                         f"ORIGIN has {len(seq)}\n")
    return len(seq), feats, seq


def cassette_span(feats):
    """Packaged AAV genome = first ITR start .. last ITR end."""
    itrs = [f for f in feats if "ITR" in f[0].upper()]
    if len(itrs) < 2:
        sys.exit("FATAL: fewer than two ITR features found")
    return min(f[1] for f in itrs), max(f[2] for f in itrs)


def load_fasta(path):
    d, cur = {}, None
    for line in open(path):
        if line.startswith(">"):
            cur = line[1:].split()[0]; d[cur] = []
        elif cur:
            d[cur].append(line.strip())
    return {k: "".join(v).upper() for k, v in d.items()}


def shared_label(cassettes, k):
    """For each construct, label every base by the set of constructs sharing the
    k-mer window that starts there (either strand). Trailing k-1 bases inherit
    the last window's label."""
    index = collections.defaultdict(set)
    for name, s in cassettes.items():
        for st in (s, rc(s)):
            for i in range(len(st) - k + 1):
                index[st[i:i + k]].add(name)
    labels = {}
    for name, s in cassettes.items():
        lab = []
        for i in range(len(s) - k + 1):
            lab.append(frozenset(index[s[i:i + k]]))
        lab += [lab[-1]] * (k - 1) if lab else []
        labels[name] = lab
    return labels, index


def blocks_from_labels(lab):
    out, start = [], 0
    for i in range(1, len(lab) + 1):
        if i == len(lab) or lab[i] != lab[start]:
            out.append((start + 1, i, lab[start]))   # 1-based inclusive
            start = i
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    cassettes, feat_map, spans, plasmid_len = {}, {}, {}, {}

    for reporter, fn in CONSTRUCTS:
        path = os.path.join(PLASMID, fn)
        L, feats, seq = parse_genbank(path)
        cs, ce = cassette_span(feats)
        cassettes[reporter] = seq[cs - 1:ce]
        spans[reporter] = (cs, ce, fn)
        plasmid_len[reporter] = L
        # shift features into cassette coordinates, keep those fully inside
        feat_map[reporter] = [(k, s - cs + 1, e - cs + 1, st)
                              for (k, s, e, st) in feats if s >= cs and e <= ce]

    # ---- cassette FASTA ----------------------------------------------------
    with open(os.path.join(OUT, "cassettes.fa"), "w") as fh:
        for r, s in cassettes.items():
            fh.write(f">{r}_cassette\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")

    # ---- features ----------------------------------------------------------
    with open(os.path.join(OUT, "features.tsv"), "w") as fh:
        fh.write("cassette\tfeature\tstart\tend\tlength\tstrand\n")
        for r, fl in feat_map.items():
            for (k, s, e, st) in sorted(fl, key=lambda x: x[1]):
                fh.write(f"{r}\t{k}\t{s}\t{e}\t{e - s + 1}\t{st}\n")

    # ---- where the v1 contigs sit inside the cassettes ----------------------
    v1 = load_fasta(V1FA)
    v1_hits = {}          # v1 contig -> list of (cassette, start, end)
    modelled = {r: [False] * len(s) for r, s in cassettes.items()}
    for cname, cseq in v1.items():
        for r, s in cassettes.items():
            p = s.find(cseq)
            while p != -1:
                v1_hits.setdefault(cname, []).append((r, p + 1, p + len(cseq)))
                for i in range(p, p + len(cseq)):
                    modelled[r][i] = True
                p = s.find(cseq, p + 1)
            # also try the reverse complement
            p = s.find(rc(cseq))
            while p != -1:
                v1_hits.setdefault(cname, []).append((r, p + 1, p + len(cseq), "rc"))
                for i in range(p, p + len(cseq)):
                    modelled[r][i] = True
                p = s.find(rc(cseq), p + 1)

    # ---- sharing blocks at both k -----------------------------------------
    lab31, idx31 = shared_label(cassettes, K_INDEX)
    lab90, _     = shared_label(cassettes, K_READ)

    with open(os.path.join(OUT, "blocks.tsv"), "w") as fh:
        fh.write("cassette\tstart\tend\tlength\tshared_with\tn_constructs\t"
                 "in_v1_reference\tfeatures\n")
        for r in cassettes:
            for (s, e, share) in blocks_from_labels(lab90[r]):
                ov = sorted({k for (k, fs, fe, st) in feat_map[r]
                             if fs <= e and fe >= s})
                mod = sum(modelled[r][s - 1:e])
                tag = "yes" if mod == (e - s + 1) else ("partial" if mod else "no")
                fh.write(f"{r}\t{s}\t{e}\t{e - s + 1}\t"
                         f"{','.join(sorted(share))}\t{len(share)}\t{tag}\t"
                         f"{','.join(ov) if ov else '-'}\n")

    # ---- k-mer index for stage B ------------------------------------------
    with open(os.path.join(OUT, "cassette_kmers.tsv"), "w") as fh:
        for kmer, names in idx31.items():
            fh.write(f"{kmer}\t{','.join(sorted(names))}\n")

    # ---- per-base mask handed to stage B -----------------------------------
    with open(os.path.join(OUT, "cassette_basemask.tsv"), "w") as fh:
        fh.write("cassette\tpos\tin_v1_reference\tn_constructs_k90\tfeature\n")
        for r in cassettes:
            fl = sorted(feat_map[r], key=lambda x: (x[2] - x[1]))
            for i in range(len(cassettes[r])):
                feat = "-"
                for (k, fs, fe, st) in fl:          # smallest containing feature wins
                    if fs <= i + 1 <= fe:
                        feat = k; break
                fh.write(f"{r}\t{i+1}\t{int(modelled[r][i])}\t{len(lab90[r][i])}\t{feat}\n")

    # ---- summary -----------------------------------------------------------
    with open(os.path.join(OUT, "summary.txt"), "w") as fh:
        def w(*a):
            print(*a); fh.write(" ".join(str(x) for x in a) + "\n")

        w("=== cassettes reconstructed from the plasmid maps ===")
        w(f"{'reporter':<10} {'plasmid':<9} {'plasmid_bp':>10} {'cassette':>12} {'cassette_bp':>12}")
        for r, _ in CONSTRUCTS:
            cs, ce, fn = spans[r]
            w(f"{r:<10} {fn:<9} {plasmid_len[r]:>10} {f'{cs}..{ce}':>12} {len(cassettes[r]):>12}")

        w("")
        w("=== how much of each cassette the v1 reference actually contains ===")
        w(f"{'reporter':<10} {'cassette_bp':>12} {'modelled':>10} {'missing':>9} {'pct_missing':>12}")
        tot_c = tot_m = 0
        for r, _ in CONSTRUCTS:
            L = len(cassettes[r]); mod = sum(modelled[r]); miss = L - mod
            tot_c += L; tot_m += mod
            w(f"{r:<10} {L:>12} {mod:>10} {miss:>9} {100.0 * miss / L:>11.1f}%")
        w(f"{'TOTAL':<10} {tot_c:>12} {tot_m:>10} {tot_c - tot_m:>9} "
          f"{100.0 * (tot_c - tot_m) / tot_c:>11.1f}%")

        w("")
        w("=== v1 contigs located in the cassettes (exact match) ===")
        for cname, cseq in sorted(v1.items()):
            hits = v1_hits.get(cname)
            if hits:
                for h in hits:
                    w(f"  {cname:<24} len {len(cseq):>5}  ->  {h[0]}_cassette {h[1]}..{h[2]}"
                      + ("  (reverse complement)" if len(h) > 3 else ""))
            else:
                w(f"  {cname:<24} len {len(cseq):>5}  ->  NOT FOUND exactly in any cassette")

        w("")
        w("=== unmodelled stretches inside each cassette (>=20 bp) ===")
        for r, _ in CONSTRUCTS:
            runs, st = [], None
            for i, m in enumerate(modelled[r]):
                if not m and st is None: st = i
                elif m and st is not None:
                    runs.append((st + 1, i)); st = None
            if st is not None: runs.append((st + 1, len(modelled[r])))
            w(f"  -- {r}_cassette")
            for (s, e) in runs:
                if e - s + 1 < 20: continue
                ov = sorted({k for (k, fs, fe, sd) in feat_map[r] if fs <= e and fe >= s})
                w(f"     {s:>6}..{e:<6} {e - s + 1:>6} bp   {','.join(ov) if ov else '(no annotated feature)'}")

        w("")
        w(f"=== cassette sequence shared between constructs (window = {K_READ} bp = R2 length) ===")
        w("A read falling in a multi-construct block cannot be assigned uniquely.")
        for r, _ in CONSTRUCTS:
            L = len(cassettes[r])
            cnt = collections.Counter(len(x) for x in lab90[r])
            w(f"  {r}_cassette ({L} bp): "
              + "  ".join(f"{n}-construct {c} bp ({100.0*c/L:.1f}%)"
                          for n, c in sorted(cnt.items())))

        w("")
        w("=== is the v1 'shared' contig really shared? ===")
        w("The v1 AAV_WPRE_SV40pA_shared contig was cut from p24. Checking whether")
        w("the equivalent region of the other two constructs is identical to it.")
        sh = v1.get("AAV_WPRE_SV40pA_shared")
        if sh:
            for r, _ in CONSTRUCTS:
                s = cassettes[r]
                if s.find(sh) != -1:
                    w(f"  {r:<10} contains the v1 shared contig exactly at {s.find(sh)+1}")
                else:
                    # longest prefix of the shared contig present in this cassette
                    lo, hi = 0, len(sh)
                    while lo < hi:
                        mid = (lo + hi + 1) // 2
                        if s.find(sh[:mid]) != -1: lo = mid
                        else: hi = mid - 1
                    w(f"  {r:<10} does NOT contain it; longest matching prefix = {lo} bp "
                      f"of {len(sh)}")

        w("")
        w(f"stage A outputs in {OUT}")

if __name__ == "__main__":
    main()
