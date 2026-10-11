#!/usr/bin/env python3
r"""EXPERIMENTAL: collocation-count scan for the unprojected, even and odd systems.

Random exact rational points above the unit hemisphere, midpoint floating least squares only.
Writes global_cusp_symmetry_scan_result.json.  Not a certificate; flips no flag.
"""
import numpy as np, random
from fractions import Fraction
from flint import arb, ctx
import global_quarter_turn as gq
from track_b_two_cusp_hejhal import ValidatedWhittaker, assemble_physical_rows_from_points, midpoint_matrix
from global_cusp_symmetry_experiment import sector_matrix, R_FROZEN
ctx.prec=128
def points(count, seed):
    rng=random.Random(seed); out=[]
    while len(out)<count:
        x1=Fraction(rng.randint(-20,20),40); x2=Fraction(rng.randint(-20,20),40); y=Fraction(rng.randint(30,70),40)
        if x1*x1+x2*x2+y*y>1: out.append({"x1":str(x1),"x2":str(x2),"y":str(y)})
    return out
def solve(M, pts, eps):
    b=ValidatedWhittaker(arb(R_FROZEN)); rows,led,mi,mz=assemble_physical_rows_from_points(M,pts,b,None)
    a=midpoint_matrix(rows); n=a.shape[1]
    B=np.eye(n) if eps is None else sector_matrix(M,eps,n)
    A=a@B
    nidx=next(k for k,m in enumerate(mi) if m[:2]==(1,0)); nc=int(np.nonzero(B[nidx])[0][0])
    rest=[j for j in range(A.shape[1]) if j!=nc]
    sol,*_=np.linalg.lstsq(A[:,rest],-A[:,nc],rcond=1e-13)
    c=np.zeros(A.shape[1],dtype=complex); c[nc]=1; c[rest]=sol
    full=B@c
    keys=[('i',m[0],m[1]) for m in mi]+[('z',m[0],m[1]) for m in mz]
    sv=np.linalg.svd(np.vstack([A,np.eye(A.shape[1])[[nc]]])/np.linalg.norm(A,axis=0).clip(1e-300),compute_uv=False)
    rank=int(np.linalg.matrix_rank(np.vstack([A,np.eye(A.shape[1])[[nc]]]),tol=max(A.shape)*np.finfo(float).eps*np.linalg.norm(A,2)))
    return dict(zip(keys,full)), float(np.linalg.norm(A@c)), rank, A.shape[1]
def main():
    import json, sys
    from pathlib import Path
    counts=[int(x) for x in sys.argv[1:]] or [6,24,60]
    result={"tier":"EXPERIMENTAL (floating lstsq on Arb midpoints)","r":R_FROZEN,"bits":128,"seed":7,"scans":[]}
    for npts in counts:
        pts=points(npts,7)
        for name,eps in (("odd",-1),("even",1),("unprojected",None)):
            result["scans"].append(run_one(npts,pts,name,eps))
    Path("global_cusp_symmetry_scan_result.json").write_text(json.dumps(result,indent=2),encoding="utf-8")
def run_one(npts,pts,name,eps):
    sols={}; info={}
    for M in (1,2,3,4,5,6):
        s,res,rank,nun=solve(M,pts,eps); sols[M]=s; info[M]=(res,rank,nun)
    drift=[]
    for M in (1,2,3,4,5):
        common=[k for k in sols[M] if k in sols[M+1]]
        drift.append(np.linalg.norm(np.array([sols[M][k]-sols[M+1][k] for k in common])))
    rec={"points":npts,"system":name,"drift_M_to_M+1_l2":[float(d) for d in drift],
         "lstsq_residual_M1_to_M6":[info[M][0] for M in range(1,7)],"rank_at_M6":info[6][1],"unknowns_at_M6":info[6][2]}
    print(npts,name,"drift",["%.2e"%d for d in rec["drift_M_to_M+1_l2"]],"rank/unk@M6 %d/%d"%(rec["rank_at_M6"],rec["unknowns_at_M6"]))
    return rec
if __name__=="__main__":
    main()
