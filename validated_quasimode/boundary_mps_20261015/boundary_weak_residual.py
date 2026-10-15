#!/usr/bin/env python3
"""Picard Q-odd boundary weak-residual experiment.

NUMERICAL DIAGNOSTICS ONLY for J and Chebyshev coefficients.
Exact rational arithmetic validates the conditional functional-analytic gates.
A rigorous interval 2D DCT certificate remains outstanding.
"""
import csv
import json
import math
from fractions import Fraction as Q
from pathlib import Path

import numpy as np
from scipy.fft import dct
from numpy.polynomial.chebyshev import chebval
from numpy.polynomial.legendre import leggauss
from sys import path as sys_path
HERE=Path(__file__).resolve().parent
sys_path.insert(0,str(HERE/'inputs'))
from level_one_quasimode_analytic import modes_and_orbits

R=Q('6.622119348586011'); LAM=1+R*R
MU_LOWER=Q(13,1000000)     # Provided by prior independently replayed MPFR slab proof
BOUNDARY_L2_GATE=Q(4,10**12)
BOUNDARY_LINF_GATE=Q(4,10**12)
EPS=Q(1,100000)


def exact_conditions():
    # For y0>=1/sqrt(2), δ=1/3, Ctrace <= 81/4.
    # Area: (1-x²-t²)^(-3/2) <=1+4(x²+t²), x²+t²<=1/2.
    # For eigenvalues, 1+LAM+EPS <49 and d>=EPS implies
    # ||r||V* >= EPS/sqrt(1+LAM+EPS)||f||.
    assert LAM<Q(45)
    assert 1+LAM+EPS <49
    assert Q(9,2)*Q(4,10**12)*Q(7,1) < MU_LOWER*EPS
    assert Q(13,1)**2/Q(24,1)**2 + Q(25,1)**2/Q(48,1)**2 < 1-Q(13,20)**2
    assert Q(144)*Q(2*10**7)*Q(2,3)**128 < Q(1,10**13)
    assert Q(13,10)*37*(Q(11,10)*Q(13,400)*3780000 + Q(13,10)*Q(17,10)*120000) < 2*10**7
    return {
      'eigenvalue_lambda':str(LAM),'trace_constant_upper':str(Q(81,4)),
      'boundary_area_upper':str(Q(5,6)),
      'conditional_L2_gate':str(BOUNDARY_L2_GATE),
      'conditional_Linfty_gate':str(BOUNDARY_LINF_GATE),
      'eigenvalue_radius':str(EPS),
      'two_dimensional_interpolation_tail_upper':str(Q(144)*Q(2*10**7)*Q(2,3)**128),
      'polyellipse_M_upper':'20000000',
      'status':'CONDITIONAL; certified nodal Chebyshev coefficient intervals MISSING'
    }


def build_basis():
    dat=json.loads((HERE/'inputs'/'analytic_quasimode_coefficients_50_reflection_even.json').read_text())
    modes,orbits,S=modes_and_orbits(50)
    mp={}
    for idx,o in enumerate(dat['orbits']):
        a=Q(o['coefficient_real'])
        for mode_id in orbits[idx]:
            ab=modes[mode_id][:2]
            mp[ab]=a*int(S[idx,mode_id])
    terms=[]
    for a,b,n in modes:
        if a>b>=0 and mp[(a,b)]!=0:
            c=(2 if b==0 else 4)*mp[(a,b)]
            terms.append((a,b,n,c))
    assert len(terms)==21
    assert sum(abs(v[-1]) for v in terms)<37
    assert max(a+b for a,b,_,_ in terms)<=9
    assert max(a for a,_,_,_ in terms)<=7
    return terms


def radial_tables():
    pol={};drv={}
    radial_path=HERE/'inputs'/'radial_coeffs_all.csv'
    if not radial_path.exists():
        radial_path=HERE.parent/'interval_taylor2_local'/'inputs'/'radial'/'radial_coeffs_all.csv'
    with open(radial_path,newline='') as fh:
        for n,k,lo,hi,dlo,dhi in csv.reader(fh):
            if n.startswith('#'):continue
            n,k=int(n),int(k)
            pol.setdefault(n,np.zeros(121))[k]=float((Q(lo)+Q(hi))/2)
            drv.setdefault(n,np.zeros(121))[k]=float((Q(dlo)+Q(dhi))/2)
    assert len(pol)==24 and all(len(v)==121 for v in pol.values())
    return pol,drv


def numerical_boundary(x,t,terms,pol,drv):
    y=np.sqrt(1-x*x-t*t)
    tau=(15/8)*(y-17/15)
    P=2*np.pi
    J=np.zeros_like(y);f=np.zeros_like(y)
    for a,b,n,c_exact in terms:
        c=float(c_exact)
        F=chebval(tau,pol[n]); Fy=chebval(tau,drv[n])
        ax=P*a*x;bx=P*b*x;at=P*a*t;bt=P*b*t
        ca,cb,ct,dt=np.cos(ax),np.cos(bx),np.cos(at),np.cos(bt)
        T=ca*dt-cb*ct
        Tx=P*(-a*np.sin(ax)*dt+b*np.sin(bx)*ct)
        Tt=P*(-b*ca*np.sin(bt)+a*cb*np.sin(at))
        f+=c*F*T
        J-=y*c*(F*(x*Tx+t*Tt)+y*Fy*T)
    return y,J,f


def direct_gauss(nx,nt,terms,pol,drv):
    xs,wx=leggauss(nx);ts,wt=leggauss(nt)
    x=.5*xs[:,None]+np.zeros((1,nt))
    t=np.zeros((nx,1))+.25*(ts[None,:]+1)
    y,J,f=numerical_boundary(x,t,terms,pol,drv)
    w=.5*wx[:,None]*.25*wt[None,:]/y**3
    return {'nx':nx,'nt':nt,'J_L2_numerical':float(np.sqrt(np.sum(w*J*J))),
      'max_sampled_abs_J':float(np.max(np.abs(J))),
      'one_sided_floor_flux_numerical':float(np.sum(w*J*f)),
      'estimated_eigenvalue_inclusion_radius_using_numeric_J':float(7*(9/2)*np.sqrt(np.sum(w*J*J))/float(MU_LOWER)),
      'status':'NONRIGOROUS numerical quadrature'}


def cheb_coefficients(n,terms,pol,drv):
    k=np.arange(n)
    a=np.cos(np.pi*(k+.5)/n)
    x=.5*a[:,None]+np.zeros((1,n))
    t=np.zeros((n,1))+.25*(1+a[None,:])
    _,J,_=numerical_boundary(x,t,terms,pol,drv)
    cc=dct(dct(J,type=2,axis=0),type=2,axis=1)/(n*n)
    cc[0,:]/=2; cc[:,0]/=2
    return {'grid':n,'max_abs_nodal_J':float(np.max(np.abs(J))),
      'abs_coefficient_sum_numerical':float(np.sum(np.abs(cc))),
      'high_order_coefficient_sum_numerical':float(np.sum(np.abs(cc[25:,:]))+np.sum(np.abs(cc[:25,25:]))),
      'interval_tail_bound_if_polyellipse_20million_proved':float(Q(144)*Q(2*10**7)*Q(2,3)**n),
      'status':'NUMERICAL DCT; not certified coefficient intervals'},cc


def run():
    terms=build_basis();pol,drv=radial_tables()
    report={'exact_conditional_theorem':exact_conditions(),
      'gauss':[direct_gauss(n,t,terms,pol,drv) for n,t in [(32,28),(60,50),(90,70),(120,96)]],
      'chebyshev':[]}
    for n in [64,96,128]:
        row,_=cheb_coefficients(n,terms,pol,drv)
        report['chebyshev'].append(row)
    print(json.dumps(report,indent=2))
    out=HERE/'reports'/'boundary_diagnostics.json';out.write_text(json.dumps(report,indent=2)+'\n')
    return report

if __name__=='__main__':run()
