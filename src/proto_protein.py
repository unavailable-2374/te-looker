#!/usr/bin/env python3
"""
Direction #5: protein-anchored ancient/divergent TE discovery. Anchors = TE-protein
(RepeatPeps) loci from diamond blastx on the unmasked residual (diverged LINE/L1 etc.
invisible to exact 16-mers). Extract element windows around each anchor -> cd-hit cluster
-> strand-consistent adaptive spoa POA consensus. Downstream copy-validates + re-masks;
RepeatMasker then masks the FULL element + all diverged copies (incl. those w/o protein).
"""
import time, subprocess, re, numpy as np
from collections import defaultdict
GENOME="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa"
REGIONS="/tmp/rp_regions.bed"
CDHIT="/home/shuoc/tool/miniconda3/envs/PGTA/bin/cd-hit-est"
SPOA="/home/shuoc/tool/miniconda3/envs/PGTA/bin/spoa"
OUT="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/proto"
CORE=400; W0=800; WMAX=8000; AGREE=0.65; MINCOV=0.50; POA_MAX=40; GAPTOL=25; MINFAM=3
B={65:0,67:1,71:2,84:3}
def rc(s): return s.translate(str.maketrans("ACGTNacgtn","TGCANtgcan"))[::-1]
def is_tandem(seq):
    a=np.frombuffer(seq.encode(),np.uint8); n=len(a)
    if n<100: return False
    return any(np.mean(a[:n-p]==a[p:])>0.8 for p in list(range(2,51))+[75,100,150,200,300,400,500] if p<n)
def poa(seqs):
    tmp=f"{OUT}/_poa_p.fa"
    with open(tmp,"w") as f:
        for i,s in enumerate(seqs): f.write(f">{i}\n{s}\n")
    msa=subprocess.run([SPOA,"-r","1",tmp],capture_output=True,text=True).stdout
    rows=[];cur=[]
    for ln in msa.splitlines():
        if ln.startswith('>'):
            if cur: rows.append("".join(cur)); cur=[]
        else: cur.append(ln.strip())
    if cur: rows.append("".join(cur))
    if len(rows)<3: return None,False
    L=len(rows[0]); arr=np.array([[ord(c) for c in r] for r in rows])
    cons=[];ag=np.zeros(L);cov=np.zeros(L)
    for j in range(L):
        col=arr[:,j]; ng=col[(col!=ord('-'))&(col!=ord('N'))]; cov[j]=len(ng)/len(rows)
        if len(ng)==0: cons.append('-'); continue
        v,c=np.unique(ng,return_counts=True); cons.append(chr(v[np.argmax(c)])); ag[j]=c.max()/len(ng)
    good=(cov>=MINCOV)&(ag>=AGREE)
    runs=[];st=None;last=None
    for j in range(L):
        if good[j]:
            if st is None: st=j
            last=j
        elif st is not None and j-last>GAPTOL: runs.append((st,last)); st=None
    if st is not None: runs.append((st,last))
    if not runs: return None,False
    lo,hi=max(runs,key=lambda r:r[1]-r[0])
    consensus="".join(c for c in (cons[j] for j in range(lo,hi+1)) if c!='-')
    sat=good[:30].mean()>0.5 or good[-30:].mean()>0.5
    return consensus,sat
def main():
    t0=time.time()
    chroms={};nm=None;buf=[]
    for ln in open(GENOME):
        if ln[0]=='>':
            if nm: chroms[nm]="".join(buf)
            nm=ln[1:].split()[0];buf=[]
        else: buf.append(ln.strip())
    chroms[nm]="".join(buf)
    loci=[]
    for ln in open(REGIONS):
        c,s,e=ln.split()[:3]; loci.append((c,(int(s)+int(e))//2))
    print(f"[anchors] {len(loci)} protein-anchored loci ({time.time()-t0:.1f}s)")
    core=[]
    for (c,mid) in loci:
        if mid-CORE<0 or mid+CORE>len(chroms[c]): continue
        w=chroms[c][mid-CORE:mid+CORE].upper()
        if w.count('N')>0.2*len(w): continue
        core.append((c,mid,w))
    cf=f"{OUT}/prot_core.fa"
    with open(cf,"w") as f:
        for i,(c,mid,w) in enumerate(core): f.write(f">{i}\n{w}\n")
    subprocess.run([CDHIT,"-i",cf,"-o",f"{OUT}/prot_clust","-c","0.80","-n","5",
                    "-aS","0.5","-r","1","-d","0","-T","64","-M","32000","-g","1"],
                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    clusters=[];cur=None
    for ln in open(f"{OUT}/prot_clust.clstr"):
        if ln[0]=='>':
            if cur: clusters.append(cur)
            cur=[]
        else:
            mid=int(re.search(r'>(\d+)\.\.\.',ln).group(1))
            cur.append((mid,'+' if '*' in ln else re.search(r'at ([+-])',ln).group(1)))
    if cur: clusters.append(cur)
    fams=[c for c in clusters if len(c)>=MINFAM]
    print(f"[cluster] {len(clusters)} clusters, {len(fams)} families>=3 ({time.time()-t0:.1f}s)")
    def extract(mem,w):
        out=[]
        for i,strnd in mem:
            c,mid,_=core[i]; a=mid-w;b=mid+w
            if a<0 or b>len(chroms[c]): continue
            s=chroms[c][a:b].upper()
            if s.count('N')>0.3*len(s): continue
            s=s if strnd=='+' else rc(s)
            out.append(s.replace('N','A'))
        return out
    outc=open(f"{OUT}/protein_consensi.fa","w"); nf=0
    for c in fams:
        mem=c[:]
        if len(mem)>POA_MAX: mem=mem[:POA_MAX]
        w=W0; cons=None
        while True:
            seqs=extract(mem,w)
            if len(seqs)<3: break
            cn,sat=poa(seqs)
            if cn is None: break
            cons=cn
            if sat and 2*w<=WMAX and len(cn)<7000 and not is_tandem(cn): w*=2; continue
            break
        if not cons or len(cons)<80: continue
        nf+=1
        outc.write(f">protein_fam_{nf} members={len(c)} win={w} len={len(cons)}\n")
        for i in range(0,len(cons),80): outc.write(cons[i:i+80]+"\n")
    outc.close()
    print(f"[done] {nf} protein-anchored consensi; wall={time.time()-t0:.1f}s")
if __name__=="__main__": main()
