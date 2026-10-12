#!/usr/bin/env python3
"""Independent Arb jet reconstruction and certified nonzero-minor rank lower bounds.

Exact geometric row equivalence supplies the universal upper bound 14.
This does not certify the numerical kernel, a trial, or an eigenvalue.
"""
import json
from pathlib import Path

from flint import arb, acb, acb_mat, ctx
import numpy as np
from scipy.linalg import qr

import ford_vertex_equivariance as geom
import global_quarter_turn as gq
from track_b_two_cusp_data import cusp_modes, mode_frequency
from track_b_two_cusp_hejhal import ValidatedWhittaker
from ford_equivariant_collocation import exact_zero_jet_columns


def reconstruct(cutoff,r):
    mi,mz=cusp_modes(cutoff); n=len(mi)+len(mz)
    backend=ValidatedWhittaker(arb(r)); jets=[]
    for x,z,t in geom.VERTICES.values():
        ax,az,at=(arb(v.numerator)/v.denominator for v in (x,z,t)); y=at.sqrt()
        j=[[acb(0) for k in range(n)] for i in range(24)]
        radial={}
        for c in range(6):
            modes=mi if c==0 else mz; off=0 if c==0 else len(mi)
            for k,mode in enumerate(modes):
                u,v,mag=mode_frequency('infinity' if c==0 else 'zero',mode)
                key=(c==0,mode[2])
                if key not in radial:
                    value,_=backend.radial(mag,y); derivative,_=backend.radial_y_derivative(mag,y)
                    radial[key]=value,derivative
                value,derivative=radial[key]
                phase=acb(0,2*arb.pi()*(u*(ax+(0 if c==0 else c-1))+v*az)).exp()
                j[c][off+k]=value*phase
                j[6+3*c][off+k]=acb(0,2*arb.pi()*u)*value*phase
                j[7+3*c][off+k]=acb(0,2*arb.pi()*v)*value*phase
                j[8+3*c][off+k]=derivative*phase/(2*y)
        jets.extend(j)
    constraints=[]
    for row in geom.reduced_vertex_constraint_rows():
        constraints.append([sum((acb(arb(v.numerator)/v.denominator)*jets[k][j]
                                for k,v in row.items()),acb(0)) for j in range(n)])
    reduced=[]
    for row in constraints:
        reduced.append([sum((row[k]*v for k,v in b.items()),acb(0))
                        for b in gq.sector_basis(cutoff,-1)])
    return reduced,backend.summary()


def rank_lower(cutoff,r):
    ctx.prec=192
    matrix,work=reconstruct(cutoff,r)
    mid=np.array([[complex(float(v.real.mid()),float(v.imag.mid())) for v in row] for row in matrix])
    zeros=exact_zero_jet_columns(cutoff)
    active=[k for k in range(mid.shape[1]) if k not in zeros]
    scaled=mid[:,active]/np.linalg.norm(mid[:,active],axis=0).clip(1e-300)
    nr=np.linalg.norm(scaled,axis=1); scaled=scaled/nr.clip(1e-300)[:,None]
    s=np.linalg.svd(scaled,compute_uv=False)
    rank=int(np.sum(s>1e-10*s[0]))
    _,_,cols=qr(scaled,pivoting=True,mode='economic')
    _,_,rows=qr(scaled[:,cols[:rank]].T,pivoting=True,mode='economic')
    rows=list(map(int,rows[:rank])); cols=[active[int(k)] for k in cols[:rank]]
    minor=acb_mat([[matrix[i][j] for j in cols] for i in rows])
    determinant=minor.det()
    lower=rank if determinant.is_finite() and not determinant.contains(0) else 0
    upper=min(14,len(matrix[0])-len(zeros))
    return {'M':cutoff,'r':r,'arb_bits':192,'selected_rows':rows,'selected_columns':cols,
            'determinant':{'real':str(determinant.real),'imag':str(determinant.imag)},
            'certified_rank_lower':lower,'exact_rank_upper':upper,'exact_rank_certified':lower==upper,
            'exact_zero_jet_columns':list(zeros),
            'numerical_rank_diagnostic':rank,'reconstructed_Arb_jet_entries':True,
            'numerical_matrix_entries_used_for_determinant':False,'work':work,
            'constraint_kernel_certified':False,'spectral_certificate':False}


if __name__=='__main__':
    results=[]
    for r in ('6.62208','6.7439020359331625'):
        for M in range(2,9):
            rec=rank_lower(M,r); results.append(rec)
            print(M,r,rec['certified_rank_lower'],rec['exact_rank_upper'],flush=True)
    out={'schema':'ford-constraint-rank-arb/v1','results':results,
         'certified':'nonzero Arb minors only; geometric rank upper14 proved separately',
         'independent_from_midpoint_jet_assembler':True,'hejhal_existence_certified':False}
    Path('ford_constraint_rank_arb.json').write_text(json.dumps(out,indent=2)+'\n')
