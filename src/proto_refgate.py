#!/usr/bin/env python3
"""
Refiner-style per-family acceptance gate. For each new-method family, compare the
ORIGINAL vs REFINED consensus by genuine genomic coverage = merged genomic bp hit at
>=80% identity over >=50% of THAT consensus's length. Accept the refined consensus only
if it covers MORE genome by this criterion. The >=50%-of-own-length rule rejects chimeras
(a bad extension dilutes the covered fraction of the longer consensus) and rewards genuine
full-length extensions. Build the gated library; downstream RM measures + copy-validates.
"""
import subprocess, re
from collections import defaultdict
GEN="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/genome/genome.fa"
BL="/home/shuoc/tool/miniconda3/envs/PGTA/bin/blastn"
LIB="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/mask_eval/lib_scaled3.fa"
REF="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/proto/refined_newmethod.fa"
OUTDIR="/scratch/shuoc/TE/Arabidopsis_thaliana/demo/mask_eval"
def read_fa(f):
    d={};n=None;s=[]
    for ln in open(f):
        if ln[0]=='>':
            if n:d[n]="".join(s)
            n=ln.split()[0][1:];s=[]
        else:s.append(ln.strip())
    if n:d[n]="".join(s)
    return d
def coverage(fa):
    """merged genomic bp per query at >=80% id, aln>=50% of query length."""
    out=subprocess.run([BL,"-query",fa,"-subject",GEN,"-evalue","1e-5","-num_threads","32",
                        "-outfmt","6 qseqid sseqid sstart send pident length qlen"],
                        capture_output=True,text=True).stdout
    iv=defaultdict(list)
    for ln in out.splitlines():
        q,s,ss,se,pid,length,qlen=ln.split()
        if float(pid)>=80 and int(length)>=0.5*int(qlen):
            a,b=sorted((int(ss),int(se))); iv[(q,s)].append((a,b))
    cov=defaultdict(int)
    for (q,s),lst in iv.items():
        lst.sort(); ce=-1
        for a,b in lst:
            if a>ce: cov[q]+=b-a+1; ce=b
            elif b>ce: cov[q]+=b-ce; ce=b
    return cov
def main():
    lib=read_fa(LIB); ref=read_fa(REF)
    nm=[n for n in ref]                          # new-method families that were refined
    # write orig subset (same names) for fair blastn
    with open(f"{OUTDIR}/_orig_nm.fa","w") as f:
        for n in nm:
            if n in lib: f.write(f">{n}\n{lib[n]}\n")
    cov_o=coverage(f"{OUTDIR}/_orig_nm.fa")
    cov_r=coverage(REF)
    accept=set(); gain=0
    for n in nm:
        if cov_r.get(n,0) > cov_o.get(n,0): accept.add(n); gain+=cov_r[n]-cov_o.get(n,0)
    print(f"families refined: {len(nm)} | refined ACCEPTED (cover more): {len(accept)} | rejected (kept original): {len(nm)-len(accept)}")
    print(f"summed coverage gain from accepted: {gain:,} bp (pre-RM proxy)")
    # build gated library: mdl-repeat + (refined if accepted else original) new-method
    with open(f"{OUTDIR}/lib_refgated.fa","w") as f:
        for n,s in lib.items():
            if n in accept: f.write(f">{n}\n{ref[n]}\n")
            else: f.write(f">{n}\n{s}\n")
    print("lib_refgated.fa written")
if __name__=="__main__": main()
