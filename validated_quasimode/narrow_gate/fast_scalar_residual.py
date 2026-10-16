#!/usr/bin/env python3
"""Scalar specialization of the 40-orbit analytic residual.
Numerical regression and integration aid; NOT certified interval arithmetic.
"""
import numpy as np
from level_one_quasimode_analytic import cutoff

class ScalarResidual:
    def __init__(self,q,orbit_coefficients):
        self.q=q
        self.b=np.asarray(q.S.T@orbit_coefficients,dtype=complex)
        self.bx=2j*np.pi*q.mu1*self.b
        self.bt=2j*np.pi*q.mu2*self.b
    def seed(self,x,t,y):
        q=self.q
        rr=np.empty((len(q.modes),len(y)))
        dr=np.empty_like(rr)
        for n, ids in q.radial_slots.items():
            val,grad=q.rad.evaluate(n,y)
            rr[ids,:]=val
            dr[ids,:]=grad
        ph=np.exp(2j*np.pi*(q.mu1[:,None]*x[None,:]+q.mu2[:,None]*t[None,:]))
        f=np.sum(self.b[:,None]*rr*ph,axis=0)
        fx=np.sum(self.bx[:,None]*rr*ph,axis=0)
        ft=np.sum(self.bt[:,None]*rr*ph,axis=0)
        fy=np.sum(self.b[:,None]*dr*ph,axis=0)
        return f,(fx,ft,fy)
    def value_and_residual(self,x,t,y):
        q=self.q
        n=len(y);G=np.zeros(n);LG=np.zeros(n);gG=np.zeros((3,n));
        N=np.zeros(n,complex);gN=np.zeros((3,n),complex)
        S=np.zeros(n,complex);T=np.zeros(n,complex)
        for gamma in q.gammas:
            xp,tp,h,zj,hj=q.transform(gamma,x,t,y)
            w,wp,wpp=cutoff(h,q.ya,q.yb)
            ids=np.flatnonzero(w>0)
            if not len(ids):continue
            f,(fx,ft,fy)=self.seed(xp[ids],tp[ids],h[ids]);
            hx,ht,hy=(v[ids] for v in hj)
            jx,jt,jy=(v[ids] for v in zj)
            gf=np.stack([fx*jx.real+ft*jx.imag+fy*hx,
                         fx*jt.real+ft*jt.imag+fy*ht,
                         fx*jy.real+ft*jy.imag+fy*hy])
            dw=np.stack([wp[ids]*hx,wp[ids]*ht,wp[ids]*hy])
            lw=h[ids]*wp[ids]-h[ids]**2*wpp[ids]
            G[ids]+=w[ids];LG[ids]+=lw;gG[:,ids]+=dw
            N[ids]+=w[ids]*f
            gN[:,ids]+=w[ids]*gf+dw*f[None,:]
            S[ids]+=lw*f
            T[ids]+=y[ids]**2*np.sum(dw*gf,axis=0)
        assert min(G)>.99999
        W=N/G
        gradW=(gN-W[None,:]*gG)/G[None,:]
        R=(S-W*LG-2*T+2*y*y*np.sum(gG*gradW,axis=0))/G
        return W,R
