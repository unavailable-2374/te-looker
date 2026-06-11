#!/usr/bin/env python3
"""
Scaled new-method discovery TARGETING mdl-repeat's gaps: run on the residual genome
(mdl-repeat-masked regions set to N) over ALL high-copy residual seeds (not 40).
Produces a library of repeats mdl-repeat MISSED -> added to mdl-repeat, raises mask ratio.
Same machinery as proto_e3e4 (seed->cd-hit cluster->adaptive-window spoa POA->boundary,
§E6 polyA preserve, §E1 runaway guard). Goal = masking coverage, so tandem/satellite KEPT.
"""
import time, subprocess, re, numpy as np
from collections import defaultdict
GENOME="/tmp/residual2.fa"; JFDB="/tmp/res16_2.jf"
JF="/home/shuoc/tool/jellyfish/bin/jellyfish"
CDHIT="/home/shuoc/tool/miniconda3/envs/PGTA/bin/cd-hit-est"
SPOA="/home/shuoc/tool/miniconda3/envs/PGTA/bin/spoa"
OUT="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/proto"
K=16; SEED_MINCOUNT=20; C_OCC=120; CAP_INST=30; MINFAM=10; CORE=300
W0=400; WMAX=8000; AGREE=0.70; MINCOV=0.50; POA_MAX=40; GAPTOL=20
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
def rc(s): return s.translate(str.maketrans("ACGTN","TGCAN"))[::-1]
def is_tandem(seq):
    a=np.frombuffer(seq.encode(),np.uint8); n=len(a)
    if n<100: return False
    for p in list(range(2,51))+[75,100,150,200,300,400,500,600]:
        if p>=n: break
        if np.mean(a[:n-p]==a[p:])>0.80: return True
    return False
def poa_consensus(seqs):
    tmp=f"{OUT}/_poa_s2.fa"
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
    L=len(rows[0]); arr=np.array([[ord(ch) for ch in r] for r in rows])
    cons=[];agree=np.zeros(L);cover=np.zeros(L)
    for j in range(L):
        col=arr[:,j]; nong=col[col!=ord('-')]; cover[j]=len(nong)/len(rows)
        if len(nong)==0: cons.append('-'); continue
        vals,cnt=np.unique(nong,return_counts=True); cons.append(chr(vals[np.argmax(cnt)])); agree[j]=cnt.max()/len(nong)
    good=(cover>=MINCOV)&(agree>=AGREE)
    runs=[];st=None;last=None
    for j in range(L):
        if good[j]:
            if st is None: st=j
            last=j
        elif st is not None and j-last>GAPTOL: runs.append((st,last)); st=None
    if st is not None: runs.append((st,last))
    if not runs: return None,False
    lo,hi=max(runs,key=lambda r:r[1]-r[0])
    while hi+1<L and cons[hi+1]=='A': hi+=1
    while lo-1>=0 and cons[lo-1]=='T': lo-=1
    consensus="".join(ch for ch in (cons[j] for j in range(lo,hi+1)) if ch!='-')
    sat=good[:30].mean()>0.5 or good[-30:].mean()>0.5
    return consensus,sat
def main():
    t0=time.time()
    chroms={};nm=None;buf=[]
    for ln in open(GENOME):
        if ln[0]=='>':
            if nm: chroms[nm]="".join(buf).upper()
            nm=ln[1:].split()[0];buf=[]
        else: buf.append(ln.strip())
    chroms[nm]="".join(buf).upper(); order=list(chroms)
    rows=[l.split() for l in subprocess.run([JF,"dump","-c","-L",str(SEED_MINCOUNT),JFDB],
            capture_output=True,text=True).stdout.splitlines() if l]
    def s2c(s):
        c=0
        for ch in s: c=(c<<2)|B[ord(ch)]
        return np.uint32(c)
    def lowcomplex(k):
        return max(k.count(b) for b in 'ACGT')>=11 or len(set(k))<=2 or any(b*8 in k for b in 'ACGT')
    seeds={s2c(r[0]):int(r[1]) for r in rows if not lowcomplex(r[0])}
    print(f"[seeds] {len(seeds)} residual seeds count>={SEED_MINCOUNT} ({time.time()-t0:.1f}s)")
    canon={};strand={};valid={}
    for c in order:
        b2=enc(chroms[c]); fwd,v=fwd_codes(b2); rc2=rc32(fwd)
        canon[c]=np.minimum(fwd,rc2); strand[c]=(fwd<=rc2); valid[c]=v
    print(f"[index] residual canonical codes built ({time.time()-t0:.1f}s)")
    sc_arr=np.array(sorted(seeds.keys()),dtype=np.uint32)
    cidx={c:i for i,c in enumerate(order)}
    per_seed=defaultdict(list)               # seed_code -> [(c,p,strand)]
    for c in order:
        cc=canon[c]; vv=valid[c]; st=strand[c]
        idx=np.nonzero(np.isin(cc,sc_arr)&vv)[0]
        for p in idx: per_seed[int(cc[p])].append((c,int(p),bool(st[p])))
    insts=[];core=[]
    for scode,occ in per_seed.items():       # §B: deterministic bottom-k cap per seed
        if len(occ)>CAP_INST:
            occ.sort(key=lambda o:((cidx[o[0]]*100000000+o[1])*2654435761)&0xffffffff)
            occ=occ[:CAP_INST]
        for (c,p,stt) in occ:
            if p-CORE<0 or p+K+CORE>len(chroms[c]): continue
            w=chroms[c][p-CORE:p+K+CORE]
            if 'N' in w: continue
            core.append(w if stt else rc(w)); insts.append((c,p,stt))
    M=len(insts); print(f"[instances] {M} from {len(per_seed)} seeds (cap {CAP_INST}/seed) ({time.time()-t0:.1f}s)")
    cf=f"{OUT}/scaled2_core.fa"
    with open(cf,"w") as f:
        for i,s in enumerate(core): f.write(f">{i}\n{s}\n")
    subprocess.run([CDHIT,"-i",cf,"-o",f"{OUT}/scaled2_clust","-c","0.80","-n","5",
                    "-aS","0.6","-r","1","-d","0","-T","64","-M","48000","-g","1"],
                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    clusters=[];cur=None
    for ln in open(f"{OUT}/scaled2_clust.clstr"):
        if ln[0]=='>':
            if cur: clusters.append(cur)
            cur=[]
        else:
            mid=int(re.search(r'>(\d+)\.\.\.',ln).group(1))
            cur.append((mid,'+' if '*' in ln else re.search(r'at ([+-])',ln).group(1)))
    if cur: clusters.append(cur)
    fams=[c for c in clusters if len(c)>=MINFAM]
    print(f"[cluster] {len(clusters)} clusters, {len(fams)} families>=3 ({time.time()-t0:.1f}s)")
    def extract(members,w):
        out=[]
        for i,strnd in members:
            c,p,st=insts[i]; a=p-w;b=p+K+w
            if a<0 or b>len(chroms[c]): continue
            s=chroms[c][a:b]
            if s.count('N')>0.3*len(s): continue
            s=s if st else rc(s)
            if strnd=='-': s=rc(s)
            out.append(s.replace('N','A'))
        return out
    outc=open(f"{OUT}/scaled2_consensi.fa","w"); nf=0
    for c in fams:
        mem=c[:]
        if len(mem)>POA_MAX:
            mem.sort(key=lambda x:x[0]); mem=mem[:POA_MAX]
        w=W0; cons=None
        while True:
            seqs=extract(mem,w)
            if len(seqs)<3: break
            cn,sat=poa_consensus(seqs)
            if cn is None: break
            cons=cn
            if sat and 2*w<=WMAX and len(cn)<6000 and not is_tandem(cn): w*=2; continue
            break
        if not cons or len(cons)<80: continue
        nf+=1
        outc.write(f">scaled2_fam_{nf} members={len(c)} win={w} len={len(cons)}\n")
        for i in range(0,len(cons),80): outc.write(cons[i:i+80]+"\n")
    outc.close()
    print(f"[done] {nf} new-method consensi (residual); wall={time.time()-t0:.1f}s")
if __name__=="__main__": main()
