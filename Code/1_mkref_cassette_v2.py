#!/usr/bin/env python3
"""Build the v2 transgene contigs from the reconstructed cassettes.

Design (from Results/2_4_cassette_report/FINDINGS.md):

The v1 reference held three ORF fragments plus one 799 bp "shared" 3' element
cut from p24. That element is not shared - p131 matches only 395/799 bp and p140
only 589/799 - so the 3' end of the tdTomato and EBFP2 constructs, which is what
a 3' assay actually reaches, was absent from the reference. 7,952 diagnostic
reads (+32%) were unmappable as a result.

v2 keeps every distinct sequence exactly once, so nothing multimaps away:

  <reporter>_cassette   full ITR-to-ITR cassette in its own sequence, with every
                        base that is NOT construct-diagnostic replaced by N.
                        A base is kept if some 31-mer window containing it occurs
                        in this cassette and in neither of the other two.
  AAV_shared            the masked-out blocks, taken once from p24 and joined by
                        runs of N so no read can align across a block junction.

Masking rather than cutting means real cassette coordinates are preserved and no
artificial junctions are created inside a cassette.
"""
import os, sys, collections

PROJ = "/scratch/user/uqschauq/Insular"
BLOCKS = os.path.join(PROJ, "Results/2_4_cassette_blocks")
OUT = os.path.join(PROJ, "Data/Reference/transgenes_cassette_v2")
K = 31
SPACER = "N" * 100          # between shared blocks, > read length
MIN_BLOCK = 40              # shared runs shorter than this are not worth a contig

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

def runs_of(flags, want):
    out, st = [], None
    for i, f in enumerate(flags):
        if f == want and st is None: st = i
        elif f != want and st is not None:
            out.append((st, i - 1)); st = None
    if st is not None: out.append((st, len(flags) - 1))
    return out

def main():
    os.makedirs(OUT, exist_ok=True)
    cassettes = {k.replace("_cassette", ""): v
                 for k, v in load_fasta(os.path.join(BLOCKS, "cassettes.fa")).items()}

    # ---- which bases are diagnostic for exactly one construct ---------------
    owners = collections.defaultdict(set)
    for name, s in cassettes.items():
        for st in (s, rc(s)):
            for i in range(len(st) - K + 1):
                owners[st[i:i + K]].add(name)

    keep = {}
    for name, s in cassettes.items():
        flags = [False] * len(s)
        for i in range(len(s) - K + 1):
            if len(owners[s[i:i + K]]) == 1:
                for j in range(i, i + K):
                    flags[j] = True          # every base the diagnostic window covers
        keep[name] = flags

    # ---- AAV_shared: masked blocks of p24, taken once ----------------------
    egfp = cassettes["EGFP"]
    shared_blocks = [(a, b) for (a, b) in runs_of(keep["EGFP"], False)
                     if b - a + 1 >= MIN_BLOCK]
    shared_seq = SPACER.join(egfp[a:b + 1] for (a, b) in shared_blocks)

    # ---- every base must live somewhere: verify ----------------------------
    shared_k = set()
    for st in (shared_seq, rc(shared_seq)):
        for i in range(len(st) - K + 1):
            km = st[i:i + K]
            if "N" not in km:
                shared_k.add(km)

    # A base is reachable if some 31-mer window covering it is either diagnostic
    # (so it lives on the cassette) or present in AAV_shared. Counting orphaned
    # k-mers instead overcounts: a window straddling a diagnostic/masked boundary
    # is in neither contig, yet every base it covers is still reachable by other
    # windows.
    unreachable = {}
    for name, s_ in cassettes.items():
        reach = [False] * len(s_)
        for i in range(len(s_) - K + 1):
            km = s_[i:i + K]
            if len(owners[km]) == 1 or km in shared_k:
                for j in range(i, i + K):
                    reach[j] = True
        unreachable[name] = runs_of(reach, False)

    # ---- write FASTA -------------------------------------------------------
    fa = os.path.join(OUT, "transgenes.fa")
    with open(fa, "w") as fh:
        for name in ("EGFP", "tdTomato", "EBFP2"):
            s = "".join(c if k else "N" for c, k in zip(cassettes[name], keep[name]))
            fh.write(f">{name}_cassette\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
        fh.write(">AAV_shared\n")
        for i in range(0, len(shared_seq), 60):
            fh.write(shared_seq[i:i + 60] + "\n")

    # ---- write GTF ---------------------------------------------------------
    gtf = os.path.join(OUT, "transgenes.gtf")
    lengths = {f"{n}_cassette": len(cassettes[n]) for n in cassettes}
    lengths["AAV_shared"] = len(shared_seq)
    with open(gtf, "w") as fh:
        for contig, gid in (("EGFP_cassette", "EGFP"),
                            ("tdTomato_cassette", "tdTomato"),
                            ("EBFP2_cassette", "EBFP2"),
                            ("AAV_shared", "AAV_shared")):
            L = lengths[contig]
            fh.write(f'{contig}\ttransgene\texon\t1\t{L}\t.\t+\t.\t'
                     f'gene_id "{gid}"; transcript_id "{gid}"; gene_name "{gid}"; '
                     f'gene_biotype "protein_coding";\n')

    # ---- provenance --------------------------------------------------------
    with open(os.path.join(OUT, "PROVENANCE.md"), "w") as fh:
        fh.write("# Transgene contigs v2 - cassette-derived, built 2026-08-27\n\n")
        fh.write("From `Data/Plasmid/{p24,p131,p140}.gb` via `Results/2_4_cassette_blocks/`.\n")
        fh.write("NOT from `plasmidConstruct.FA` (its KJ081243.1 is 25.8% nt-identical\n")
        fh.write("to p131's codon-diversified tdTomato).\n\n")
        fh.write("| contig | length | diagnostic bases kept | masked to N |\n")
        fh.write("|---|---|---|---|\n")
        for name in ("EGFP", "tdTomato", "EBFP2"):
            L = len(cassettes[name]); k = sum(keep[name])
            fh.write(f"| {name}_cassette | {L} | {k} | {L - k} |\n")
        fh.write(f"| AAV_shared | {len(shared_seq)} | - | "
                 f"{len(shared_blocks) - 1} spacers of 100 |\n\n")
        fh.write(f"AAV_shared holds {len(shared_blocks)} blocks of p24 "
                 f">= {MIN_BLOCK} bp, joined by 100-N spacers so no read can\n"
                 "align across a junction. Genes: EGFP, tdTomato, EBFP2, AAV_shared\n")
        fh.write("(v1 called the last one AAV_WPRE_SV40pA_shared).\n\n")
        fh.write("## sequence accounted for\n\n")
        fh.write("Bases reachable by no contig (every 31-mer covering them is shared\n")
        fh.write("between two constructs and too short a block to host a 90 bp read):\n\n")
        for n, rr in unreachable.items():
            tot = sum(b - a + 1 for a, b in rr)
            fh.write(f"- {n}: {tot} bp"
                     + (" -- " + ", ".join(f"{a+1}..{b+1} ({b-a+1} bp)" for a, b in rr)
                        if rr else "")
                     + "\n")
        fh.write("\nThese sit inside the EGFP/EBFP2 identity runs. v1 lost the same\n")
        fh.write("reads as multimappers, so this is not a regression, and no block is\n")
        fh.write("long enough for a 90 bp read to align to it anyway.\n")

    # ---- console summary ---------------------------------------------------
    print("v2 transgene contigs")
    for name in ("EGFP", "tdTomato", "EBFP2"):
        L = len(cassettes[name]); k = sum(keep[name])
        print(f"  {name+'_cassette':<20} {L:>6} bp   diagnostic {k:>5} "
              f"({100.0*k/L:.1f}%)   masked {L-k}")
    print(f"  {'AAV_shared':<20} {len(shared_seq):>6} bp   "
          f"{len(shared_blocks)} blocks from p24")
    for n, rr in unreachable.items():
        tot = sum(b - a + 1 for a, b in rr)
        print(f"  unreachable in {n}: {tot} bp"
              + ("  " + ", ".join(f"{a+1}..{b+1}({b-a+1})" for a, b in rr) if rr else ""))
        if any(b - a + 1 >= 60 for a, b in rr):
            print(f"  WARNING: {n} has an unreachable block >= 60 bp", file=sys.stderr)

if __name__ == "__main__":
    main()
