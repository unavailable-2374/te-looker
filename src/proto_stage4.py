#!/usr/bin/env python3
"""
Stage-4 (cross-seed clustering) + POA increment over proto_a1.
 same k=16 seeds -> oriented instances -> cd-hit-est clustering on tight seed-core
 windows (merges same-family instances across seeds, splits heterogeneous occurrences
 of one seed; clustering on the core avoids the random-flank problem) -> per-family
 strand-consistent spoa POA of the FULL windows + flank-agreement boundary trim.
Controlled comparison: SAME 40 seeds as proto_a1 (baseline 5/40 full-length).
NOTE: minimap2 base-level all-vs-all was tried first and was pathologically slow on
repeat cores (anchor explosion) -> switched to cd-hit greedy clustering.
"""
import time, subprocess, os, re, numpy as np
from collections import defaultdict
GENOME="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa"
JF="/home/shuoc/tool/jellyfish/bin/jellyfish"; JFDB="/tmp/at_k16.jf"
CDHIT="/home/shuoc/tool/miniconda3/envs/PGTA/bin/cd-hit-est"
SPOA="/home/shuoc/tool/miniconda3/envs/PGTA/bin/spoa"
OUT="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/proto"
K=16; C_OCC=200; CAP_INST=80; FLANK=1500; CORE=300; N_SEEDS=40
AGREE=0.70; MINCOV=0.50; POA_MAX=50
B={65:0,67:1,71:2,84:3}
def enc(s):
    a=np.frombuffer(s.encode(),np.uint8); c=np.full(a.shape,255,np.uint8)
    for ch,v in B.items(): c[a==ch]=v
    return c
def rc32(codes):
    out=np.zeros_like(codes); c=codes.copy()
    for _ in range(K): out=(out<<np.uint32(2))|((np.uint32(3)-(c&np.uint32(3)))&np.uint32(3)); c>>=np.uint32(2)
    return out
def fwd_codes(b2):
    n=b2.shape[0]; valid=(b2!=255); bb=b2.astype(np.uint32); bb[~valid]=0
    code=np.zeros(n,np.uint32)
    for j in range(K): code[:n-K+1]=(code[:n-K+1]<<np.uint32(2))|bb[j:j+(n-K+1)]
    code[n-K+1:]=0
    cv=np.cumsum(valid.astype(np.int64)); ok=np.zeros(n,bool)
    ok[:n-K+1]=(cv[K-1:]-np.concatenate(([0],cv[:n-K])))==K
    return code,ok
def revcomp(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def main():
    t0=time.time()
    chroms={}; nm=None; buf=[]
    for ln in open(GENOME):
        if ln[0]=='>':
            if nm: chroms[nm]="".join(buf).upper()
            nm=ln[1:].split()[0]; buf=[]
        else: buf.append(ln.strip())
    chroms[nm]="".join(buf).upper(); order=list(chroms)
    rows=[l.split() for l in subprocess.run([JF,"dump","-c","-L",str(C_OCC),JFDB],
            capture_output=True,text=True).stdout.splitlines() if l]
    rows.sort(key=lambda r:-int(r[1])); idx=np.linspace(0,len(rows)-1,N_SEEDS).astype(int)
    def s2c(s):
        c=0
        for ch in s: c=(c<<2)|B[ord(ch)]
        return np.uint32(c)
    seeds={s2c(rows[i][0]):i for i in idx}
    canon={}; strand={}; valid={}
    for c in order:
        b2=enc(chroms[c]); fwd,v=fwd_codes(b2); rc=rc32(fwd)
        canon[c]=np.minimum(fwd,rc); strand[c]=(fwd<=rc); valid[c]=v
    full=[]; core=[]; seed_of=[]
    for sc,si in seeds.items():
        occ=[]
        for c in order:
            for p in np.nonzero((canon[c]==sc)&valid[c])[0]: occ.append((c,int(p),bool(strand[c][p])))
        if len(occ)>CAP_INST:
            occ.sort(key=lambda o:hash((o[0],o[1]))&0xffffffffffffffff); occ=occ[:CAP_INST]
        for (c,p,st) in occ:
            a=p-FLANK; b=p+K+FLANK
            if a<0 or b>len(chroms[c]): continue
            w=chroms[c][a:b]; w=w if st else revcomp(w)
            full.append(w); core.append(w[FLANK-CORE:FLANK+K+CORE]); seed_of.append(si)
    M=len(full); print(f"[instances] {M} oriented windows ({time.time()-t0:.1f}s)")
    cf=f"{OUT}/inst_core.fa"
    with open(cf,"w") as f:
        for i,s in enumerate(core): f.write(f">{i}\n{s}\n")
    # Stage-4: cd-hit-est greedy clustering on cores
    subprocess.run([CDHIT,"-i",cf,"-o",f"{OUT}/core_clust","-c","0.80","-n","5",
                    "-aS","0.6","-r","1","-d","0","-T","64","-M","32000","-g","1"],
                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    # parse .clstr  ->  clusters of (id, strand)
    clusters=[]; cur=None
    for ln in open(f"{OUT}/core_clust.clstr"):
        if ln[0]=='>':
            if cur: clusters.append(cur)
            cur=[]
        else:
            mid=int(re.search(r'>(\d+)\.\.\.',ln).group(1))
            if '*' in ln: cur.append((mid,'+'))
            else: cur.append((mid, re.search(r'at ([+-])',ln).group(1)))
    if cur: clusters.append(cur)
    fams=[c for c in clusters if len(c)>=3]
    seeds_per_fam=[len(set(seed_of[i] for i,_ in c)) for c in fams]
    cof={}
    for ci,c in enumerate(clusters):
        for i,_ in c: cof[i]=ci
    fps=defaultdict(set)
    for i in range(M): fps[seed_of[i]].add(cof[i])
    print(f"[stage4] {len(clusters)} clusters, {len(fams)} families(>=3); "
          f"merged(seeds/fam>1): {sum(1 for s in seeds_per_fam if s>1)}; "
          f"split(seed in>1 fam): {sum(1 for s in fps.values() if len(s)>1)} ({time.time()-t0:.1f}s)")
    # Stage-5: per family POA (full windows, reoriented by cd-hit strand) + boundary trim
    outc=open(f"{OUT}/stage4_consensi.fa","w"); nf=0
    for c in fams:
        mem=c[:]
        if len(mem)>POA_MAX:
            mem.sort(key=lambda x:hash(x[0])&0xffff); mem=mem[:POA_MAX]
        tmp=f"{OUT}/_poa.fa"
        with open(tmp,"w") as f:
            for i,strnd in mem:
                s=full[i]; s=s if strnd=='+' else revcomp(s)
                f.write(f">{i}\n{s}\n")
        msa=subprocess.run([SPOA,"-r","1",tmp],capture_output=True,text=True).stdout
        rows2=[]; curseq=[]
        for ln in msa.splitlines():
            if ln.startswith('>'):
                if curseq: rows2.append("".join(curseq)); curseq=[]
            else: curseq.append(ln.strip())
        if curseq: rows2.append("".join(curseq))
        if len(rows2)<3: continue
        L=len(rows2[0]); arr=np.array([[ord(ch) for ch in r] for r in rows2])
        cons=[]; agree=np.zeros(L); cover=np.zeros(L)
        for j in range(L):
            col=arr[:,j]; nong=col[col!=ord('-')]; cover[j]=len(nong)/len(rows2)
            if len(nong)==0: cons.append('-'); continue
            vals,cnt=np.unique(nong,return_counts=True); cons.append(chr(vals[np.argmax(cnt)]))
            agree[j]=cnt.max()/len(nong)
        good=(cover>=MINCOV)&(agree>=AGREE)
        # largest CONTIGUOUS well-supported block (tolerate short gaps) = the element,
        # excluding coincidentally-good columns out in the divergent flanks
        GAPTOL=20; runs=[]; st=None; last=None
        for j in range(L):
            if good[j]:
                if st is None: st=j
                last=j
            elif st is not None and j-last>GAPTOL:
                runs.append((st,last)); st=None
        if st is not None: runs.append((st,last))
        if not runs: continue
        lo,hi=max(runs,key=lambda r:r[1]-r[0])
        consensus="".join(ch for ch in (cons[j] for j in range(lo,hi+1)) if ch!='-')
        if len(consensus)<80: continue
        nf+=1
        outc.write(f">s4_fam_{nf} members={len(c)} poa={len(mem)} len={len(consensus)}\n")
        for i in range(0,len(consensus),80): outc.write(consensus[i:i+80]+"\n")
    outc.close()
    print(f"[done] {nf} consensi; wall={time.time()-t0:.1f}s")
if __name__=="__main__": main()
