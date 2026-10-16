#!/usr/bin/env python3
"""Independent *numerical* defect-only commutator for Picard finite quasimodes.

Never interprets Chebyshev Bessel approximation or sampled grids as interval
proof. The identity itself is exact for the underlying Fourier-Whittaker seed.
"""
import json
import math
import numpy as np
from pathlib import Path
import sys
# Run either from the flat replay archive or from the nested GitHub directory.
sys.path.insert(0,str(Path(__file__).resolve().parents[2]))
from level_one_quasimode_analytic import Quasimode, cutoff
from level_one_quasimode_verify import load_coeffs


def defects(q, x, t, y, coeff):
    """Return (W, residual) without N,G quotient differentiation or ∇W.

    Set d_g=f(gP)-f(P), H=sum w_g*d_g and W=f+H/G. Each d_g is
    an exact (L-lambda)-solution. For positive L=-y²Δ+y∂y:

      (L-lambda)W = (Σ[(Lw)d-2⟨∇w,∇d⟩])/G
          -H*LG/G² -2H|∇G|²/G³+2⟨∇H,∇G⟩/G².
    """
    n = len(x)
    f0, gf0 = q.f_and_grad(x, t, y)
    f0 = coeff @ f0
    df0 = np.asarray([coeff @ v for v in gf0])
    H = np.zeros(n, dtype=complex)
    dH = np.zeros((3,n), dtype=complex)
    G = np.ones(n)
    dG = np.zeros((3,n))
    LG = np.zeros(n)
    LHminus = np.zeros(n,dtype=complex)
    for gamma in q.gammas[1:]:
        xp,tp,h,zj,hj=q.transform(gamma,x,t,y)
        w,wp,wpp=cutoff(h,q.ya,q.yb)
        idx=np.flatnonzero(w>0)
        if len(idx)==0:continue
        fg,gg=q.f_and_grad(xp[idx],tp[idx],h[idx])
        fgam=coeff@fg
        fx,ft,fy=[coeff@v for v in gg]
        jx,jt,jy=[v[idx] for v in zj]
        hx,ht,hy=[v[idx] for v in hj]
        grad=np.array([fx*jx.real+ft*jx.imag+fy*hx,
                       fx*jt.real+ft*jt.imag+fy*ht,
                       fx*jy.real+ft*jy.imag+fy*hy])
        dg=fgam-f0[idx]
        gd=grad-df0[:,idx]
        gw=np.array([wp[idx]*hx,wp[idx]*ht,wp[idx]*hy])
        lw=h[idx]*wp[idx]-h[idx]**2*wpp[idx]
        G[idx]+=w[idx]
        dG[:,idx]+=gw
        LG[idx]+=lw
        H[idx]+=w[idx]*dg
        dH[:,idx]+=gw*dg+gd*w[idx]
        LHminus[idx]+=lw*dg-2*y[idx]**2*np.sum(gw*gd,axis=0)
    dg2=y**2*np.sum(dG*dG,axis=0)
    hdG=y**2*np.sum(dH*dG,axis=0)
    R=LHminus/G-H*LG/G**2-2*H*dg2/G**3+2*hdG/G**2
    W=f0+H/G
    return W,R


def run(coefpath,degree=118):
    src=json.loads(Path(coefpath).read_text())
    q=Quasimode(M=src['M'],r=float(src['r']),degree=degree)
    co=load_coeffs(q,coefpath)
    rng=np.random.default_rng(1167)
    x=rng.uniform(-.49,.49,42)
    t=rng.uniform(.001,.499,42)
    floor=np.sqrt(1-x*x-t*t)
    y=floor+rng.uniform(.005,.6,42)
    W,R=defects(q,x,t,y,co)
    Wcols,Rcols,_=q.W_residual_cols(x,t,y)
    a=Wcols@co;b=Rcols@co
    diff_W=float(np.max(np.abs(W-a)))
    diff_R=float(np.max(np.abs(R-b)))
    maxR=float(np.max(np.abs(R)))
    print(json.dumps({'points':len(x),'degree':degree,'max_abs_W_difference':diff_W,
         'max_abs_residual_difference':diff_R,'max_sampled_abs_residual':maxR,
         'r':q.r,'lambda':q.lam,'nonidentity_cosets':len(q.gammas)-1},indent=2))
    assert diff_R < 5e-16, (diff_R,maxR)
    assert diff_W < 5e-16, diff_W
    return {'points':len(x),'max_abs_W_difference':diff_W,
            'max_abs_residual_difference':diff_R,'max_sampled_abs_residual':maxR}

if __name__=='__main__':
    import sys
    run(sys.argv[1] if len(sys.argv)>1 else str((Path(__file__).resolve().parent/'analytic_quasimode_coefficients_50.json') if (Path(__file__).resolve().parent/'analytic_quasimode_coefficients_50.json').exists() else Path(__file__).resolve().parents[2]/'analytic_quasimode_coefficients_50.json'))
