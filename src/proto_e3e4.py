#!/usr/bin/env python3
"""
v4 §E4 (adaptive window + LTR detection) + §E3 (short-element/SINE channel) increment,
building on proto_stage4. Per-family the POA window GROWS from 400 bp until the
cross-copy support boundary becomes internal (not saturated at the window edge) ->
short elements (SINE/MITE) resolve to their true short length, long LTR elements
(truncated by the old fixed +-1500 window) extend to full length. Then classify by
structure: polyA tail (SINE), terminal direct repeat (LTR).
Controlled comparison vs proto_stage4 (fixed +-1500): same 40 seeds.
"""
import time, subprocess, re, numpy as np
from collections import defaultdict
GENOME="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa"
JF="/home/shuoc/tool/jellyfish/bin/jellyfish"; JFDB="/tmp/at_k16.jf"
CDHIT="/home/shuoc/tool/miniconda3/envs/PGTA/bin/cd-hit-est"
SPOA="/home/shuoc/tool/miniconda3/envs/PGTA/bin/spoa"
OUT="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/proto"
K=16; C_OCC=200; CAP_INST=80; CORE=300; N_SEEDS=40
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
    """short/medium-period self-similarity -> tandem array, not a dispersed TE (§E1)."""
    a=np.frombuffer(seq.encode(),np.uint8); n=len(a)
    if n<100: return False
    for p in list(range(2,51))+[75,100,150,200,300]:
        if p>=n: break
        m=np.mean(a[:n-p]==a[p:])
        if m>0.80: return True
    return False

def poa_consensus(seqs):
    """spoa MSA -> (consensus_of_largest_supported_block, edge_saturated?)"""
    tmp=f"{OUT}/_poa.fa"
    with open(tmp,"w") as f:
        for i,s in enumerate(seqs): f.write(f">{i}\n{s}\n")
    msa=subprocess.run([SPOA,"-r","1",tmp],capture_output=True,text=True).stdout
    rows=[]; cur=[]
    for ln in msa.splitlines():
        if ln.startswith('>'):
            if cur: rows.append("".join(cur)); cur=[]
        else: cur.append(ln.strip())
    if cur: rows.append("".join(cur))
    if len(rows)<3: return None,False
    L=len(rows[0]); arr=np.array([[ord(ch) for ch in r] for r in rows])
    cons=[]; agree=np.zeros(L); cover=np.zeros(L)
    for j in range(L):
        col=arr[:,j]; nong=col[col!=ord('-')]; cover[j]=len(nong)/len(rows)
        if len(nong)==0: cons.append('-'); continue
        vals,cnt=np.unique(nong,return_counts=True); cons.append(chr(vals[np.argmax(cnt)]))
        agree[j]=cnt.max()/len(nong)
    good=(cover>=MINCOV)&(agree>=AGREE)
    runs=[]; st=None; last=None
    for j in range(L):
        if good[j]:
            if st is None: st=j
            last=j
        elif st is not None and j-last>GAPTOL: runs.append((st,last)); st=None
    if st is not None: runs.append((st,last))
    if not runs: return None,False
    lo,hi=max(runs,key=lambda r:r[1]-r[0])
    # §E6: keep low-info homopolymer terminal (polyA 3' / polyT 5') past the support
    # boundary, so SINE/LINE tails survive for §E3 classification (the exact feature
    # a coverage/agreement trim would otherwise erode)
    while hi+1<L and cons[hi+1]=='A': hi+=1
    while lo-1>=0 and cons[lo-1]=='T': lo-=1
    consensus="".join(ch for ch in (cons[j] for j in range(lo,hi+1)) if ch!='-')
    # saturated if the supported block reaches near either MSA end (element fills window)
    sat = good[:30].mean()>0.5 or good[-30:].mean()>0.5
    return consensus, sat

def classify(seq):
    n=len(seq); tags=[]
    if n>9000: tags.append("suspect_long/tandem?")   # runaway-growth guard (§E1 seg-dup/tandem)
    # SINE-ish: short + polyA/polyT run near an end (tail now preserved by §E6)
    end=seq[-60:]; start=seq[:60]
    polya = bool(re.search(r'A{6,}',end) or re.search(r'T{6,}',start))
    if n<500 and polya: tags.append("SINE?")
    elif n<500: tags.append("short")
    # LTR: terminal direct repeat (5' ~ 3', same orientation)
    ltr=0
    if n>=600:
        L0=min(500,n//3); a=seq[:L0]; b=seq[-L0:]
        best=0
        for off in range(-30,31):  # small offset tolerance, ungapped
            m=t=0
            for i in range(L0):
                j=i+off
                if 0<=j<L0:
                    t+=1; m+=(a[i]==b[j])
            if t>=100 and m/t>best: best=m/t; ltrlen=t
        if best>=0.80: ltr=ltrlen; tags.append(f"LTR({ltrlen}bp/{best:.0%})")
    return ",".join(tags) if tags else "interspersed"

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
        b2=enc(chroms[c]); fwd,v=fwd_codes(b2); rc2=rc32(fwd)
        canon[c]=np.minimum(fwd,rc2); strand[c]=(fwd<=rc2); valid[c]=v
    insts=[]; core=[]
    for sc,si in seeds.items():
        occ=[]
        for c in order:
            for p in np.nonzero((canon[c]==sc)&valid[c])[0]: occ.append((c,int(p),bool(strand[c][p])))
        if len(occ)>CAP_INST:
            occ.sort(key=lambda o:hash((o[0],o[1]))&0xffffffffffffffff); occ=occ[:CAP_INST]
        for (c,p,st) in occ:
            if p-CORE<0 or p+K+CORE>len(chroms[c]): continue
            w=chroms[c][p-CORE:p+K+CORE]; core.append(w if st else rc(w)); insts.append((c,p,st))
    M=len(insts); print(f"[instances] {M} ({time.time()-t0:.1f}s)")
    cf=f"{OUT}/inst_core.fa"
    with open(cf,"w") as f:
        for i,s in enumerate(core): f.write(f">{i}\n{s}\n")
    subprocess.run([CDHIT,"-i",cf,"-o",f"{OUT}/core_clust","-c","0.80","-n","5",
                    "-aS","0.5","-r","1","-d","0","-T","64","-M","32000","-g","1"],
                   stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
    clusters=[]; cur=None
    for ln in open(f"{OUT}/core_clust.clstr"):
        if ln[0]=='>':
            if cur: clusters.append(cur)
            cur=[]
        else:
            mid=int(re.search(r'>(\d+)\.\.\.',ln).group(1))
            cur.append((mid,'+' if '*' in ln else re.search(r'at ([+-])',ln).group(1)))
    if cur: clusters.append(cur)
    fams=[c for c in clusters if len(c)>=3]
    print(f"[stage4] {len(clusters)} clusters, {len(fams)} families ({time.time()-t0:.1f}s)")

    def extract(members,w):
        out=[]
        for i,strnd in members:
            c,p,st=insts[i]; a=p-w; b=p+K+w
            if a<0 or b>len(chroms[c]): continue
            s=chroms[c][a:b]; s=s if st else rc(s)
            if strnd=='-': s=rc(s)
            out.append(s)
        return out

    outc=open(f"{OUT}/e3e4_consensi.fa","w"); nf=0; wsizes=[]; tagcount=defaultdict(int)
    for c in fams:
        mem=c[:]
        if len(mem)>POA_MAX:
            mem.sort(key=lambda x:hash(x[0])&0xffff); mem=mem[:POA_MAX]
        w=W0; cons=None; grew=0
        while True:                                  # §E4 adaptive growth (guarded)
            seqs=extract(mem,w)
            if len(seqs)<3: break
            cn,sat=poa_consensus(seqs)
            if cn is None: break
            # §E1 runaway guard: if a doubling barely shrank edge support AND the
            # consensus is already long, the members are a tandem/seg-dup block, not a
            # dispersed TE -> stop growing (prevents the 12.8 kb / 51 GB blowup).
            tandem = cn is not None and is_tandem(cn)
            cons=cn
            if sat and 2*w<=WMAX and len(cn)<6000 and not tandem:
                w*=2; grew+=1; continue
            break
        if not cons or len(cons)<80: continue
        nf+=1; wsizes.append(w)
        tag=classify(cons); tagcount[tag.split(',')[0].split('(')[0]]+=1
        outc.write(f">e_fam_{nf} members={len(c)} win={w} len={len(cons)} struct={tag}\n")
        for i in range(0,len(cons),80): outc.write(cons[i:i+80]+"\n")
    outc.close()
    print(f"[E4] adaptive windows used: {sorted(set(wsizes))}  median={int(np.median(wsizes))}")
    print(f"[E3] structural tags: {dict(tagcount)}")
    print(f"[done] {nf} consensi; wall={time.time()-t0:.1f}s")
if __name__=="__main__": main()
