#!/usr/bin/env python3
"""Independent verification of the FROZEN level-one coefficient list.

No optimization is run here. Analytic residual is compared against a separate
finite-difference definition at hand-picked nonsingular points, and independent
Gauss/Sobol quadrature. Experimental doubles only; no interval certificate.
"""
import argparse,json,math,hashlib,time
from pathlib import Path
import numpy as np
from scipy.stats import qmc
from level_one_quasimode_analytic import Quasimode,make_grid,matrix_on_grid,measure

def load_coeffs(q,path):
    src=json.loads(Path(path).read_text());assert src['M']==q.M
    assert len(src['orbits'])==q.n_orb
    for j,d in enumerate(src['orbits']):
        assert [tuple(q.modes[i][:2]) for i in q.orbits[j]]==[tuple(x) for x in d['modes']]
    out=np.array([complex(float(c['coefficient_real']),float(c['coefficient_imag'])) for c in src['orbits']])
    assert abs(out[q.fix_idx]*q.fix_sign-1)<1e-11
    return out

def direct_fd(q,coeff,p,h):
    x,t,y=p
    def W(x,t,y):
        U,_,_=q.W_residual_cols(np.array([x]),np.array([t]),np.array([y]))
        return U[0]@coeff
    U=W(x,t,y)
    wxx=(W(x+h,t,y)-2*U+W(x-h,t,y))/h**2
    wtt=(W(x,t+h,y)-2*U+W(x,t-h,y))/h**2
    wyy=(W(x,t,y+h)-2*U+W(x,t,y-h))/h**2
    wy=(W(x,t,y+h)-W(x,t,y-h))/(2*h)
    return -y*y*(wxx+wtt+wyy)+y*wy-q.lam*U

def qmc_grid(n,seed):
    u=qmc.Sobol(3,scramble=True,seed=seed).random_base2(round(math.log2(n)))
    x=u[:,0]-.5;t=.5*u[:,1];flo=np.sqrt(1-x*x-t*t)
    y=flo+(5/3-flo)*u[:,2]
    w=.5*(5/3-flo)/(y**3)/len(x)
    return x,t,y,w

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--coefficients',required=True);ap.add_argument('--output',default='analytic_quasimode_independent_check.json')
    ap.add_argument('--degree',type=int,default=118);ap.add_argument('--gauss',default='18,14,30');ap.add_argument('--sobol',type=int,default=8192)
    args=ap.parse_args();t0=time.monotonic()
    frozen=json.loads(Path(args.coefficients).read_text())
    q=Quasimode(M=frozen['M'],r=float(frozen['r']),degree=args.degree)
    co=load_coeffs(q,args.coefficients)
    g=tuple(int(x) for x in args.gauss.split(','))
    gg=make_grid(*g);W,R,w=matrix_on_grid(q,gg)
    gmetric=measure(q,co,W,R,w,with_tail=False)
    sm=[]
    for seed in (17,):
        grid=qmc_grid(args.sobol,seed)
        U,V,weights=matrix_on_grid(q,grid)
        sm.append({'seed':seed,**measure(q,co,U,V,weights)})
    fd=[]
    points=[(.125,.16,1.04),(.22,.29,.93)]
    for p in points:
        x,t,y=p
        u,r,G=q.W_residual_cols(np.array([x]),np.array([t]),np.array([y]))
        exact=r[0]@co
        vals=[]
        for h in (.004,.002):
            v=direct_fd(q,co,p,h)
            vals.append({'h':h,'fd_real':v.real,'fd_imag':v.imag,'abs_difference':abs(v-exact)})
        fd.append({'point':p,'analytic_residual_real':exact.real,'analytic_residual_imag':exact.imag,'finite_differences':vals})
    # isometric height gradient identity test |nabla h|_H^2=h^2
    identities=[]
    x=np.array([.17,-.31]);t=np.array([.19,.41]);y=np.array([.9,1.18]);
    for gamma in q.gammas[1:]:
        _x,_t,h,_z,gradh=q.transform(gamma,x,t,y)
        ratio=(y*y*sum(v*v for v in gradh))/(h*h)
        identities.append(max(abs(ratio-1)))
    out={'tier':'EXPERIMENTAL independent analytic and finite-difference check, no certified bounds',
      'coefficients_sha256':hashlib.sha256(Path(args.coefficients).read_bytes()).hexdigest(),
      'gauss_nodes':len(w),'gauss':gmetric,'sobol_nodes':args.sobol,'sobol':sm,'fd_check':fd,
      'height_isometry_max_error':float(max(identities)),'elapsed_seconds':time.monotonic()-t0}
    Path(args.output).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('fd_check',)},indent=2),flush=True)
    print('finite difference errors:',[(p['point'],[round(v['abs_difference'],13) for v in p['finite_differences']]) for p in fd],flush=True)
if __name__=='__main__':main()
