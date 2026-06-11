#!/usr/bin/env python3
"""
Boundary completion (Refiner-style) for new-method families: rebuild each family's
consensus from its ACTUAL genomic copies (RepeatMasker .out of lib_scaled3) extracted
with generous flanks, then extend the boundary outward only while cross-copy agreement
holds (stops at unique flanks -> integrity-safe). Fuller consensi mask the full element
length + diverged copies. Output refined library; merge with mdl-repeat; re-mask.
"""
import time, subprocess, re, numpy as np
from collections import defaultdict
GENOME="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa"
RMOUT="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/mask_eval/scaled3/genome.fa.out"
LIB="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/mask_eval/lib_scaled3.fa"
SPOA="/home/shuoc/tool/miniconda3/envs/PGTA/bin/spoa"
OUT="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/proto"
FL=600; SAMPLE=30; MINCOPY=5; AGREE=0.70; MINCOV=0.50; GAPTOL=20
def rc(s): return s.translate(str.maketrans("ACGTNacgtn","TGCANtgcan"))[::-1]
def poa_refine(seqs):
    tmp=f"{OUT}/_poa_r.fa"
    with open(tmp,"w") as f:
        for i,s in enumerate(seqs): f.write(f">{i}\n{s}\n")
    msa=subprocess.run([SPOA,"-r","1",tmp],capture_output=True,text=True).stdout
    rows=[];cur=[]
    for ln in msa.splitlines():
        if ln.startswith('>'):
            if cur: rows.append("".join(cur)); cur=[]
        else: cur.append(ln.strip())
    if cur: rows.append("".join(cur))
    if len(rows)<3: return None
    L=len(rows[0]); arr=np.array([[ord(c) for c in r] for r in rows])
    cons=[];agree=np.zeros(L);cover=np.zeros(L)
    for j in range(L):
        col=arr[:,j]; nong=col[(col!=ord('-'))&(col!=ord('N'))]; cover[j]=len(nong)/len(rows)
        if len(nong)==0: cons.append('-'); continue
        v,c=np.unique(nong,return_counts=True); cons.append(chr(v[np.argmax(c)])); agree[j]=c.max()/len(nong)
    good=(cover>=MINCOV)&(agree>=AGREE)
    runs=[];st=None;last=None
    for j in range(L):
        if good[j]:
            if st is None: st=j
            last=j
        elif st is not None and j-last>GAPTOL: runs.append((st,last)); st=None
    if st is not None: runs.append((st,last))
    if not runs: return None
    lo,hi=max(runs,key=lambda r:r[1]-r[0])
    while hi+1<L and cons[hi+1]=='A': hi+=1
    while lo-1>=0 and cons[lo-1]=='T': lo-=1
    return "".join(c for c in (cons[j] for j in range(lo,hi+1)) if c!='-')
def main():
    t0=time.time()
    chroms={};nm=None;buf=[]
    for ln in open(GENOME):
        if ln[0]=='>':
            if nm: chroms[nm]="".join(buf)
            nm=ln[1:].split()[0];buf=[]
        else: buf.append(ln.strip())
    chroms[nm]="".join(buf)
    # parse RM .out -> copies per new-method family
    copies=defaultdict(list)
    for ln in open(RMOUT):
        x=ln.split()
        if len(x)<11 or not x[0].isdigit(): continue
        name=x[9]
        if not (name.startswith("scaled_fam") or name.startswith("scaled2_fam")): continue
        c=x[4]; s=int(x[5])-1; e=int(x[6]); strand='+' if x[8]=='+' else '-'
        copies[name].append((c,s,e,strand))
    print(f"[copies] {len(copies)} new-method families with RM copies ({time.time()-t0:.1f}s)")
    # load original consensi (fallback)
    orig={};n=None;s=[]
    for ln in open(LIB):
        if ln[0]=='>':
            if n: orig[n]="".join(s)
            n=ln.split()[0][1:];s=[]
        else: s.append(ln.strip())
    if n: orig[n]="".join(s)
    out=open(f"{OUT}/refined_newmethod.fa","w"); refined=0; longer=0
    for name,cps in copies.items():
        if len(cps)<MINCOPY:
            if name in orig: out.write(f">{name}\n{orig[name]}\n")
            continue
        cps.sort(key=lambda t:(t[0],t[1]))
        if len(cps)>SAMPLE:
            step=len(cps)//SAMPLE; cps=cps[::step][:SAMPLE]
        seqs=[]
        for (c,s,e,strand) in cps:
            a=max(0,s-FL); b=min(len(chroms[c]),e+FL)
            w=chroms[c][a:b].upper()
            if w.count('N')>0.3*len(w): continue
            seqs.append(w if strand=='+' else rc(w))
        if len(seqs)<3:
            if name in orig: out.write(f">{name}\n{orig[name]}\n")
            continue
        ref=poa_refine(seqs)
        base=orig.get(name,"")
        if ref and len(ref)>=max(80,len(base)):     # keep refined only if >= original (fuller)
            out.write(f">{name}\n{ref}\n"); refined+=1; longer+= (len(ref)>len(base)+20)
        elif base:
            out.write(f">{name}\n{base}\n")
    out.close()
    print(f"[refine] refined {refined} families ({longer} got longer); wall={time.time()-t0:.1f}s")
if __name__=="__main__": main()
