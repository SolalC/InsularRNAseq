#!/usr/bin/env python3
"""Step 2.4 stage C: read-level construct assignment on the reconstructed cassettes.

A read is DIAGNOSTIC for construct C if it contains at least one 31-mer that
occurs in C's cassette and in neither of the other two. This replaces the
per-base uniqueness label used in the first pass, which tagged a base by the
90 bp window starting on it and so credited bases upstream of a divergence
point as construct-specific when a read overlapping only them is not.

Each diagnostic read is then split by whether its diagnostic 31-mer already
exists in the v1 contigs:
  captured_today  - v1 could in principle assign this read
  recoverable     - the sequence is absent from v1, so the read is lost now

Reads SAM on stdin.  Usage: 2_4_cassette_diagnostic.py <sample> <stream> <outdir>
"""
import os, sys, collections

PROJ = "/scratch/user/uqschauq/Insular"
BLOCKS = os.path.join(PROJ, "Results/2_4_cassette_blocks")
V1FA = os.path.join(PROJ, "Data/Reference/transgenes/transgenes.fa")
K = 31

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

def main():
    sample, stream, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    os.makedirs(outdir, exist_ok=True)

    cassettes = {k.replace("_cassette", ""): v
                 for k, v in load_fasta(os.path.join(BLOCKS, "cassettes.fa")).items()}
    v1 = load_fasta(V1FA)

    owners = collections.defaultdict(set)
    for name, s in cassettes.items():
        for st in (s, rc(s)):
            for i in range(len(st) - K + 1):
                owners[st[i:i + K]].add(name)
    diag = {km: next(iter(v)) for km, v in owners.items() if len(v) == 1}

    v1k = set()
    for s in v1.values():
        for st in (s, rc(s)):
            for i in range(len(st) - K + 1):
                v1k.add(st[i:i + K])

    sys.stderr.write(f"{sample}/{stream}: {len(owners)} cassette 31-mers, "
                     f"{len(diag)} diagnostic for exactly one construct\n")

    n = 0
    hits = collections.Counter()
    for line in sys.stdin:
        if line[0] == "@":
            continue
        f = line.split("\t", 11)
        if len(f) < 10:
            continue
        n += 1
        s = f[9]
        found = {}
        for i in range(len(s) - K + 1):
            km = s[i:i + K]
            c = diag.get(km)
            if c:
                # True once any diagnostic k-mer for c is already in v1
                found[c] = found.get(c, False) or (km in v1k)
        if len(found) == 1:
            c, in_v1 = next(iter(found.items()))
            hits[f"diagnostic_{c}"] += 1
            hits[f"diagnostic_{c}_" + ("captured_today" if in_v1 else "recoverable")] += 1
        elif len(found) > 1:
            hits["ambiguous_multi_construct"] += 1

    out = os.path.join(outdir, f"{sample}_{stream}_diagnostic.tsv")
    with open(out, "w") as fh:
        fh.write("sample\tstream\tmetric\tvalue\n")
        fh.write(f"{sample}\t{stream}\treads_scanned\t{n}\n")
        for c in sorted(cassettes):
            for suff in ("", "_captured_today", "_recoverable"):
                key = f"diagnostic_{c}{suff}"
                fh.write(f"{sample}\t{stream}\t{key}\t{hits[key]}\n")
        fh.write(f"{sample}\t{stream}\tambiguous_multi_construct\t"
                 f"{hits['ambiguous_multi_construct']}\n")
    sys.stderr.write(f"{sample}/{stream}: scanned {n}\n")

if __name__ == "__main__":
    main()
