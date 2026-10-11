import sys, collections
K=31
FA="/scratch/user/uqschauq/Insular/Data/Reference/GRCm39_transgenes/fasta/genome.fa"
TG=("EGFP","tdTomato","EBFP2","AAV_WPRE_SV40pA_shared")
seq={}; cur=None
for line in open(FA):
    if line[0]=='>':
        cur=line[1:].split()[0]; 
        if cur in TG: seq[cur]=[]
        else: cur=None
    elif cur: seq[cur].append(line.strip())
seq={k:''.join(v).upper() for k,v in seq.items()}
comp=str.maketrans("ACGTN","TGCAN")
rc=lambda s:s.translate(comp)[::-1]
km={}
for name,s in seq.items():
    for st in (s, rc(s)):
        for i in range(len(st)-K+1):
            km.setdefault(st[i:i+K], set()).add(name)
# drop k-mers shared between contigs -> keep only diagnostic ones, but track both
uniq={k:next(iter(v)) for k,v in km.items() if len(v)==1}
sys.stderr.write(f"kmers total={len(km)} contig-unique={len(uniq)}\n")
hits=collections.Counter(); n=0
for line in sys.stdin:
    f=line.split('\t',10)
    if len(f)<10: continue
    n+=1; s=f[9]
    found=set()
    for i in range(0, len(s)-K+1):
        g=uniq.get(s[i:i+K])
        if g: found.add(g)
    for g in found: hits[g]+=1
    if found: hits['ANY']+=1
sys.stderr.write(f"reads scanned={n}\n")
print(f"reads_scanned\t{n}")
for g in TG+("ANY",): print(f"{g}\t{hits[g]}")
