#!/usr/bin/env python3
"""
Minimal prototype of v4 Stage 1-2 + A1 (extend-to-growing-consensus) + Stage-5
flank-agreement boundary, on a real genome. HONEST SCOPE:
  - exact k=16 canonical seeds (weight>=16 per V4 §A4), C_occ cap via hash bottom-k
  - per-seed extend-to-consensus = each occurrence aligned ONCE (ungapped, anchored
    at the seed) to the growing column profile  -> demonstrates A1's O(occ) cost
  - boundary by column agreement decaying to background
NOT implemented (these are the full build): Stage-4 cross-seed clustering, gapped
POA/indels, acceptance stats, classification, protein channel, nesting.
So this validates the A1 CORE MECHANISM + resource cost on real data, not a finished caller.
"""
import sys, time, subprocess, numpy as np

GENOME = "/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa"
JF     = "/home/shuoc/tool/jellyfish/bin/jellyfish"
JFDB   = "/tmp/at_k16.jf"
K      = 16
C_OCC  = 200
FLANK  = 1500          # element half-window (elements up to ~3 kb; longer LTRs truncated - known limit)
N_SEEDS= 40            # number of high-count seeds to grow families from
AGREE  = 0.70          # flank-agreement boundary threshold
MINCOV = 0.50          # min fraction of members covering a column to call it

B = {65:0,67:1,71:2,84:3}  # A C G T
def encode_chrom(seq):
    a = np.frombuffer(seq.encode(), dtype=np.uint8)
    code = np.full(a.shape, 255, np.uint8)
    for ch,v in B.items(): code[a==ch]=v
    return code  # 0..3, 255=N/other

def rc_code32(codes):
    # reverse-complement a uint32 array of 16x2bit codes (b0 high .. b15 low)
    out = np.zeros_like(codes)
    c = codes.copy()
    for _ in range(K):
        out = (out<<np.uint32(2)) | ((np.uint32(3)-(c & np.uint32(3))) & np.uint32(3))
        c >>= np.uint32(2)
    return out

def fwd_codes(b2):
    # rolling 16-mer forward codes; positions where any base is N -> invalid
    n = b2.shape[0]
    valid = (b2!=255)
    bb = b2.astype(np.uint32); bb[~valid]=0
    code = np.zeros(n, np.uint32)
    for j in range(K):
        code[:n-K+1] = (code[:n-K+1]<<np.uint32(2)) | bb[j:j+(n-K+1)]
    code[n-K+1:] = 0
    # validity of the 16-mer starting at i = all 16 bases valid
    vv = valid.astype(np.uint8)
    win = np.ones(n, np.uint8)
    cv = np.cumsum(vv)
    okcount = np.empty(n, np.int64); okcount[:n-K+1] = cv[K-1:] - np.concatenate(([0],cv[:n-K]))
    valid16 = np.zeros(n, bool); valid16[:n-K+1] = okcount[:n-K+1]==K
    return code, valid16

def revcomp_seq(s):
    t = str.maketrans("ACGTN","TGCAN"); return s.translate(t)[::-1]

def main():
    t0=time.time()
    # ---- load genome ----
    chroms={}; name=None; buf=[]
    for line in open(GENOME):
        if line[0]=='>':
            if name: chroms[name]="".join(buf).upper()
            name=line[1:].split()[0]; buf=[]
        else: buf.append(line.strip())
    if name: chroms[name]="".join(buf).upper()
    order=list(chroms);
    print(f"[load] {len(order)} chroms, {sum(len(s) for s in chroms.values())/1e6:.1f} Mb  ({time.time()-t0:.1f}s)")

    # ---- Stage 1: high-count seeds from jellyfish (count>=C_OCC) ----
    dump = subprocess.run([JF,"dump","-c","-L",str(C_OCC),JFDB],capture_output=True,text=True).stdout
    rows=[ln.split() for ln in dump.splitlines() if ln]
    rows.sort(key=lambda r:-int(r[1]))
    # sample N_SEEDS spread across the high-count list (diverse families, not just the very top)
    idx=np.linspace(0,len(rows)-1,N_SEEDS).astype(int)
    seeds=[rows[i][0] for i in idx]
    def s2code(s):
        c=0
        for ch in s: c=(c<<2)|B[ord(ch)]
        return np.uint32(c)
    seed_codes={ s2code(s):s for s in seeds }
    print(f"[stage1] {len(rows)} seeds count>={C_OCC}; grow from {len(seed_codes)} sampled  ({time.time()-t0:.1f}s)")

    # ---- per-chrom canonical codes ----
    canon={}; strand={}; valid={}
    for nm in order:
        b2=encode_chrom(chroms[nm])
        fwd,v16=fwd_codes(b2)
        rc=rc_code32(fwd)
        can=np.minimum(fwd,rc)
        st=(fwd<=rc)  # True -> '+' orientation
        canon[nm]=can; strand[nm]=st; valid[nm]=v16
    print(f"[index] canonical 16-mer codes built  ({time.time()-t0:.1f}s)")

    # ---- Stage 2 + A1: per seed, gather occ (cap via hash bottom-k), extend-to-consensus ----
    targets=np.array(sorted(seed_codes.keys()),dtype=np.uint32)
    cons_out=open("/scratch/shuoc/TE/Arabidopsis_thaliana/demo/proto/proto_consensi.fa","w")
    bed=open("/scratch/shuoc/TE/Arabidopsis_thaliana/demo/proto/proto_members.bed","w")
    nfam=0; total_align=0
    for sc,seedseq in seed_codes.items():
        occ=[]  # (chrom,pos,strand)
        for nm in order:
            hit=np.nonzero((canon[nm]==sc)&valid[nm])[0]
            for p in hit: occ.append((nm,int(p),bool(strand[nm][p])))
        if len(occ)<3: continue
        # C_occ cap via hash bottom-k (deterministic)
        if len(occ)>C_OCC:
            occ.sort(key=lambda o:hash((o[0],o[1]))& 0xffffffffffffffff)
            occ=occ[:C_OCC]
        # build oriented windows anchored at seed (seed lands at offset FLANK)
        W=2*FLANK+K; cols=np.full((len(occ),W),255,np.uint8); used=0
        for (nm,p,st) in occ:
            s=chroms[nm]; a=p-FLANK; b=p+K+FLANK
            if a<0 or b>len(s): continue
            w=s[a:b]
            if not st: w=revcomp_seq(w)
            arr=np.frombuffer(w.encode(),np.uint8)
            row=np.full(W,255,np.uint8)
            for ch,v in B.items(): row[arr==ch]=v
            cols[used]=row; used+=1
        if used<3: continue
        cols=cols[:used]; total_align+=used   # A1: each member aligned ONCE -> O(occ)
        # column majority + agreement
        cons=[]; agree=[]; cover=[]
        for j in range(W):
            col=cols[:,j]; nonN=col[col!=255]
            if len(nonN)==0: cons.append('N'); agree.append(0); cover.append(0); continue
            cnt=np.bincount(nonN,minlength=4); mb=int(np.argmax(cnt))
            cons.append("ACGT"[mb]); agree.append(cnt[mb]/len(nonN)); cover.append(len(nonN)/used)
        agree=np.array(agree); cover=np.array(cover)
        # boundary: extend from anchor (FLANK) until agreement decays below AGREE (sustained)
        def edge(rng):
            below=0; last=FLANK
            for j in rng:
                ok=(cover[j]>=MINCOV and agree[j]>=AGREE)
                if ok: last=j; below=0
                else:
                    below+=1
                    if below>=30: break
            return last
        L=edge(range(FLANK,-1,-1)); R=edge(range(FLANK,W))
        consensus="".join(cons[L:R+1]).strip("N")
        if len(consensus)<80: continue
        nfam+=1
        cons_out.write(f">proto_fam_{nfam} seed={seedseq} members={used} len={len(consensus)}\n")
        for i in range(0,len(consensus),80): cons_out.write(consensus[i:i+80]+"\n")
        for (nm,p,st) in occ:
            bed.write(f"{nm}\t{max(0,p-FLANK)}\t{p+K+FLANK}\tproto_fam_{nfam}\t.\t{'+' if st else '-'}\n")
    cons_out.close(); bed.close()
    print(f"[A1] {nfam} draft families; total member-to-consensus alignments={total_align} "
          f"(O(occ), not O(occ^2))  ({time.time()-t0:.1f}s)")
    print(f"[done] wall={time.time()-t0:.1f}s")

if __name__=="__main__":
    main()
