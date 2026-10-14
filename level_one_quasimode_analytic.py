#!/usr/bin/env python3
"""Analytic commutator check for the finite Picard odd Poincare quasimode.

No finite differences, no global mesh, no hidden numerical coefficient vector.
Numerical diagnostics ONLY: mpmath Bessel nodes -> Chebyshev approximants.
The generated rational coefficient JSON is reproducible, NOT interval certified.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
import time
from pathlib import Path
import numpy as np
from numpy.polynomial import chebyshev as C
import mpmath as mp

RSTAR = 6.622119348586011


def modes_and_orbits(M):
    k=math.isqrt(M)
    modes=sorted([(a,b,a*a+b*b) for a in range(-k,k+1) for b in range(-k,k+1) if 0<a*a+b*b<=M],key=lambda v:(v[2],v[0],v[1]))
    lookup={(a,b):i for i,(a,b,n) in enumerate(modes)}
    seen=set(); orbits=[]
    for i,(a,b,n) in enumerate(modes):
        if i in seen: continue
        orb=[lookup[(a,b)],lookup[(-b,a)],lookup[(-a,-b)],lookup[(b,-a)]]
        seen.update(orb); orbits.append(orb)
    S=np.zeros((len(orbits),len(modes)),dtype=float)
    for j,orb in enumerate(orbits): S[j,orb]=[1,-1,1,-1]
    return modes,orbits,S


def representatives():
    result=[(1+0j,0j,0j,1+0j)]
    for a in range(-3,4):
        for b in range(-3,4):
            d=complex(a,b)
            if a*a+b*b<=2:result.append((0j,-1+0j,1+0j,d))
    c=1+1j
    for d,inv in ((1,1),(-1,-1),(1j,-1j),(-1j,1j)):
        result.append((complex(inv),0j,c,complex(d)))
    return result


def cutoff(h,ya,yb):
    # C2 quintic. Endpoint derivatives must be exactly 0, without crosses.
    t=(h-ya)/(yb-ya)
    t=np.clip(t,0,1)
    q=t**3*(10-15*t+6*t*t)
    active=(h>ya)&(h<yb)
    qp=np.where(active,30*t*t*(1-t)*(1-t)/(yb-ya),0.)
    qpp=np.where(active,60*t*(1-t)*(1-2*t)/(yb-ya)**2,0.)
    return q,qp,qpp


class Radial:
    def __init__(self,r,norms,degree=120,ymin=.59,ymax=4.2,dps=24):
        self.r=r;self.names=sorted(set(int(n) for n in norms));self.ymin=ymin;self.ymax=ymax
        mp.mp.dps=dps
        kk=np.arange(degree+1)
        tt=np.cos(np.pi*(kk+.5)/(degree+1))
        ys=.5*(ymin+ymax)+.5*(ymax-ymin)*tt
        self.coef={};self.deriv={}
        t0=time.monotonic()
        for n in self.names:
            vals=np.array([float(mp.re(mp.mpf(float(y))*mp.besselk(1j*mp.mpf(str(r)),2*mp.pi*mp.sqrt(n)*mp.mpf(float(y))))) for y in ys])
            # nodes coincide with true Chebyshev-Gauss points; fit is stable.
            coef=C.chebfit(tt,vals,degree)
            self.coef[n]=coef
            self.deriv[n]=C.chebder(coef)*2/(ymax-ymin)
        self.build_seconds=time.monotonic()-t0
    def evaluate(self,n,y):
        t=(2*y-(self.ymin+self.ymax))/(self.ymax-self.ymin)
        return C.chebval(t,self.coef[n]), C.chebval(t,self.deriv[n])
    def midpoint_error_check(self):
        # Held-out radial nodes, independent of Chebyshev fitting nodes
        ys=(.601,.7,.88,1.2,1.6666666,2.1,3.41)
        results=[]
        for n in self.names:
            for y in ys:
                target=complex(mp.mpf(y)*mp.besselk(1j*mp.mpf(str(self.r)),2*mp.pi*mp.sqrt(n)*mp.mpf(y)))
                v,_=self.evaluate(n,np.array([y])); results.append(abs(v[0]-target.real))
        return max(results)


class Quasimode:
    def __init__(self,M=30,r=RSTAR,ya=.6,yb=.7,degree=120,dps=24):
        self.M,self.r,self.lam,self.ya,self.yb=M,r,1+r*r,ya,yb
        self.modes,self.orbits,self.S=modes_and_orbits(M)
        self.mu1=np.array([m[0] for m in self.modes],dtype=float)
        self.mu2=np.array([m[1] for m in self.modes],dtype=float)
        self.norms=np.array([m[2] for m in self.modes],dtype=int)
        self.rad=Radial(r,self.norms,degree=degree,dps=dps)
        self.gammas=representatives()
        self.n_orb=len(self.orbits)
        self.fix_idx=next(k for k,orb in enumerate(self.orbits) if any(self.modes[i][:2]==(1,0) for i in orb))
        self.fix_sign=self.S[self.fix_idx,next(i for i,m in enumerate(self.modes) if m[:2]==(1,0))]
        self.radial_slots={n:np.where(self.norms==n)[0] for n in self.rad.names}

    @staticmethod
    def transform(g,x,t,y):
        a,b,c,d=g
        z=x+1j*t
        q=c*z+d
        D=(q.real*q.real+q.imag*q.imag)+abs(c)**2*y*y
        p=a*z+b
        num=p*np.conj(q)+a*np.conj(c)*y*y
        zp=num/D
        h=y/D
        # Jacobian for the full upper half space action.
        Dx=2*np.real(c*np.conj(q)); Dt=2*np.real(1j*c*np.conj(q)); Dy=2*abs(c)**2*y
        nx=a*np.conj(q)+p*np.conj(c)
        nt=1j*a*np.conj(q)-1j*p*np.conj(c)
        ny=2*a*np.conj(c)*y
        zx=(nx-zp*Dx)/D;zt=(nt-zp*Dt)/D;zy=(ny-zp*Dy)/D
        hx=-h*Dx/D; ht=-h*Dt/D; hy=(1-h*Dy)/D
        return zp.real,zp.imag,h,(zx,zt,zy),(hx,ht,hy)

    def f_and_grad(self,x,t,y):
        npts=len(x)
        radial=np.empty((len(self.modes),npts)); rad_y=np.empty_like(radial)
        for n,ids in self.radial_slots.items():
            v,d=self.rad.evaluate(n,y)
            radial[ids]=v;rad_y[ids]=d
        theta=2*np.pi*(self.mu1[:,None]*x[None,:]+self.mu2[:,None]*t[None,:])
        phase=np.exp(1j*theta)
        val=self.S@(radial*phase)
        fx=self.S@((2j*np.pi*self.mu1[:,None])*radial*phase)
        ft=self.S@((2j*np.pi*self.mu2[:,None])*radial*phase)
        fy=self.S@(rad_y*phase)
        return val,(fx,ft,fy)

    def W_residual_cols(self,x,t,y):
        n=len(x);norb=self.n_orb
        G=np.zeros(n); LG=np.zeros(n); gG=np.zeros((3,n))
        N=np.zeros((norb,n),dtype=complex)
        gN=np.zeros((3,norb,n),dtype=complex)
        L_N_psi=np.zeros((norb,n),dtype=complex)
        gradpsi_gradf=np.zeros((norb,n),dtype=complex)
        for gamma in self.gammas:
            xp,tp,h,zj,hj=self.transform(gamma,x,t,y)
            w,wp,wpp=cutoff(h,self.ya,self.yb)
            active=w>0
            if not np.any(active): continue
            idx=np.flatnonzero(active)
            vals,gv=self.f_and_grad(xp[idx],tp[idx],h[idx])
            hxa,hta,hya=(v[idx] for v in hj)
            # physical 3-component gradient of f after gamma, chain rule
            jx,jt,jy=(v[idx] for v in zj)
            fx,ft,fy=gv
            gradf=np.stack([
                fx*jx.real[None,:]+ft*jx.imag[None,:]+fy*hxa[None,:],
                fx*jt.real[None,:]+ft*jt.imag[None,:]+fy*hta[None,:],
                fx*jy.real[None,:]+ft*jy.imag[None,:]+fy*hya[None,:],
            ])
            dw=np.stack([wp[idx]*hxa,wp[idx]*hta,wp[idx]*hya])
            lw=h[idx]*wp[idx]-h[idx]*h[idx]*wpp[idx]
            Wt=w[idx]
            G[idx]+=Wt
            LG[idx]+=lw
            gG[:,idx]+=dw
            N[:,idx]+=Wt[None,:]*vals
            gN[:,:,idx]+=Wt[None,None,:]*gradf+dw[:,None,:]*vals[None,:,:]
            L_N_psi[:,idx]+=lw[None,:]*vals
            gradpsi_gradf[:,idx]+=y[idx][None,:]**2*np.sum(dw[:,None,:]*gradf,axis=0)
        if np.any(G<.999999999999):
            raise ArithmeticError(f'G < 1: {float(np.min(G))}')
        W=N/G[None,:]
        gradW=(gN-W[None,:,:]*gG[:,None,:])/G[None,None,:]
        gG_gradW=np.sum(gG[:,None,:]*gradW,axis=0)*y[None,:]**2
        resid=(L_N_psi-W*LG[None,:]-2*gradpsi_gradf+2*gG_gradW)/G[None,:]
        return W.T.copy(),resid.T.copy(),G


def make_grid(nx=12,nt=10,ny=18,vertical_splits=1, x_offset=0.,t_offset=0.):
    # Independent positive tensor-product Gauss rules; y respects curved floor exactly.
    ax,wx=np.polynomial.legendre.leggauss(nx)
    at,wt=np.polynomial.legendre.leggauss(nt)
    ay,wy=np.polynomial.legendre.leggauss(ny)
    x=.5*ax; wx=.5*wx
    t=.25*(at+1);wt=.25*wt
    xx=[];tt=[];yy=[];weights=[]
    for i,xi in enumerate(x):
        for j,tj in enumerate(t):
            ymin=math.sqrt(1-xi*xi-tj*tj)
            ymax=5/3
            for k in range(vertical_splits):
                lo=ymin+(ymax-ymin)*k/vertical_splits
                hi=ymin+(ymax-ymin)*(k+1)/vertical_splits
                y=.5*(hi-lo)*ay+.5*(lo+hi)
                xx.extend([xi]*ny);tt.extend([tj]*ny);yy.extend(y)
                weights.extend(wx[i]*wt[j]*.5*(hi-lo)*wy/y**3)
    return np.array(xx),np.array(tt),np.array(yy),np.array(weights)


def matrix_on_grid(q,grid,chunk=250):
    x,t,y,weights=grid
    W=np.zeros((len(x),q.n_orb),complex);R=np.zeros_like(W)
    for s in range(0,len(x),chunk):
        W[s:s+chunk],R[s:s+chunk],G=q.W_residual_cols(x[s:s+chunk],t[s:s+chunk],y[s:s+chunk])
    return W,R,weights


def fit(q,W,R,weights,round_digits=10):
    # Analytic residual least-squares; NOT relying on almost-cancelling finite-difference stencils.
    sw=np.sqrt(weights)
    A=sw[:,None]*R
    k0=q.fix_idx;a0=q.fix_sign
    ids=[k for k in range(q.n_orb) if k!=k0]
    norms=np.linalg.norm(A[:,ids],axis=0)
    norms=np.where(norms>0,norms,1)
    coef,*_=np.linalg.lstsq(A[:,ids]/norms[None,:],-A[:,k0]*a0,rcond=1e-11)
    sol=np.zeros(q.n_orb,dtype=complex);sol[k0]=a0;sol[ids]=coef/norms
    # freeze to exact decimal rationals; retain a(1,0) = 1
    frozen=np.round(sol.real,round_digits)+1j*np.round(sol.imag,round_digits)
    return sol,frozen


def parseval_tail(q,a,Y=5/3):
    # Rigorous formula for norm on y>Y, evaluated NUMERICALLY using mpmath.
    # Area of half-square = 1/2, four Fourier modes/orbit.
    mp.mp.dps=30
    val=mp.mpf(0)
    for j,orb in enumerate(q.orbits):
        n=q.modes[orb[0]][2]
        g=lambda h: (mp.besselk(1j*mp.mpf(str(q.r)),2*mp.pi*mp.sqrt(n)*h))**2/h
        term=2*float(mp.re(mp.quad(g,[mp.mpf(str(Y)),mp.inf])))
        val+=float(abs(a[j])**2)*term
    return float(val)


def measure(q,a,W,R,weights,with_tail=False):
    Rsq=float(np.sum(weights*np.abs(R@a)**2))
    Msq=float(np.sum(weights*np.abs(W@a)**2))
    tail=parseval_tail(q,a) if with_tail else 0
    return dict(R=math.sqrt(max(Rsq,0)),mu_core=math.sqrt(max(Msq,0)),mu=math.sqrt(max(Msq+tail,0)),tail_mu=math.sqrt(tail),ratio=math.sqrt(max(Rsq,0)/max(Msq+tail,1e-300)))


def save_json(path,obj):
    text=json.dumps(obj,indent=2,sort_keys=True)+'\n'
    Path(path).write_text(text)
    return hashlib.sha256(text.encode()).hexdigest()


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--M',type=int,default=30)
    ap.add_argument('--degree',type=int,default=120)
    ap.add_argument('--nx',type=int,default=12)
    ap.add_argument('--nt',type=int,default=10)
    ap.add_argument('--ny',type=int,default=18)
    ap.add_argument('--verify-grid',default='16,14,24')
    ap.add_argument('--output',default='analytic_quasimode_experiment.json')
    ap.add_argument('--coefficients',default='analytic_quasimode_coefficients.json')
    args=ap.parse_args()
    t0=time.monotonic()
    q=Quasimode(args.M,degree=args.degree)
    print(f'Radial tables {q.rad.build_seconds:.1f}s ({len(q.rad.names)} norms), interpolation held-out absolute error ~{q.rad.midpoint_error_check():.3g}',flush=True)
    grid=make_grid(args.nx,args.nt,args.ny)
    W,R,weights=matrix_on_grid(q,grid)
    raw,frozen=fit(q,W,R,weights)
    vr=list(map(int,args.verify_grid.split(',')))
    verify=make_grid(*vr)
    V,U,vweights=matrix_on_grid(q,verify)
    tail=parseval_tail(q,frozen)
    metrics={}
    for key,a in [('fitted',raw),('decimal_frozen',frozen)]:
        metrics[key]={'train':measure(q,a,W,R,weights),'verify':measure(q,a,V,U,vweights)}
    mt=metrics['decimal_frozen']['verify']
    mt['tail_mu']=math.sqrt(tail)
    mt['mu_with_tail']=math.sqrt(mt['mu']**2+tail)
    mt['ratio_with_tail']=mt['R']/mt['mu_with_tail']
    coeffs={'tier':'EXPERIMENTAL, mpmath Chebyshev radial and analytic differential, decimal rational coefficients','r':str(q.r),'lambda':str(q.lam),'M':q.M,'ya':q.ya,'yb':q.yb,'normalization':'a(1,0) = 1','orbits':[]}
    for j,orbit in enumerate(q.orbits):
        coeffs['orbits'].append({'modes':[list(q.modes[k][:2]) for k in orbit], 'coefficient_real':format(frozen[j].real,'.10f'),'coefficient_imag':format(frozen[j].imag,'.10f')})
    sha=save_json(args.coefficients,coeffs)
    report={'tier':'EXPERIMENTAL DOUBLE PRECISION, not certified','r':q.r,'lambda':q.lam,'M':q.M,'n_orbits':q.n_orb,'cosets':len(q.gammas),'radial_degree':args.degree,'radial_heldout_max_abs':q.rad.midpoint_error_check(),'train_nodes':len(weights),'verify_nodes':len(vweights),'coefficients_sha256':sha,'metrics':metrics,'elapsed_seconds':time.monotonic()-t0}
    save_json(args.output,report)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
