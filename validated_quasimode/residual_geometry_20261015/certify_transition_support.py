#!/usr/bin/env python3
"""Exact integer/rational outer cover of compact support of Picard quasimode residual.

All inequalities use int64 with a proved small integer bound and no transcendental
functions. The result is a RIGOROUS geometric bound, NOT a residual certificate.
"""
from fractions import Fraction as Q
from pathlib import Path
import argparse,json,hashlib
import numpy as np
from level_one_quasimode_analytic import representatives
from ford_active_cosets import ACTIVE_INDICES

def distance_range(lo,hi):
    mn=np.where((lo<=0)&(hi>=0),0,np.minimum(abs(lo),abs(hi)))
    mx=np.maximum(abs(lo),abs(hi))
    return mn*mn,mx*mx

def calculate(resolution=128, y_den=192,triangle=False):
    # Common integer denominator so all comparisons are EXACT integer arithmetic.
    # Must be divisible by horizontal resolution and y denominator.
    from math import lcm
    scale=lcm(resolution,y_den)
    assert scale < 2000
    xs=np.arange(-resolution//2,resolution//2,dtype=np.int64)*(scale//resolution)
    ts=np.arange(0,resolution//2,dtype=np.int64)*(scale//resolution)
    xlo,tlo=np.meshgrid(xs,ts,indexing='ij');xlo=xlo.reshape(-1);tlo=tlo.reshape(-1)
    step=scale//resolution
    if triangle:
        inside=(xlo>=0)&(tlo<=xlo)
        xlo=xlo[inside];tlo=tlo[inside]
    xhi=xlo+step;thi=tlo+step
    xmin=np.minimum(abs(xlo),abs(xhi));xmin=np.where((xlo<=0)&(xhi>=0),0,xmin)
    xmax=np.maximum(abs(xlo),abs(xhi));tmax=np.maximum(abs(tlo),abs(thi))
    # Reject wholly outside domain using outer bounding radius.
    u_max=xmax*xmax+tmax*tmax
    n=len(xlo)
    cs=[]
    for index in ((0,1,2,4,5,13) if triangle else ACTIVE_INDICES):
        if index==0:continue
        _,_,c,d=representatives()[index];cr=int(c.real);ci=int(c.imag);dr=int(d.real);di=int(d.imag)
        # Each affine real / imaginary component varies over an interval.
        re_corners=np.stack([cr*x-ci*t+dr*scale for x in [xlo,xhi] for t in [tlo,thi]])
        im_corners=np.stack([ci*x+cr*t+di*scale for x in [xlo,xhi] for t in [tlo,thi]])
        qr_min,qr_max=distance_range(re_corners.min(axis=0),re_corners.max(axis=0))
        qi_min,qi_max=distance_range(im_corners.min(axis=0),im_corners.max(axis=0))
        cs.append((qr_min+qi_min,qr_max+qi_max,cr*cr+ci*ci,index))
    low_y=2*y_den//3;hi_y=5*y_den//3
    assert 3*low_y==2*y_den and 3*hi_y==5*y_den
    step_y=scale//y_den
    cover_counts=[];total_possible=0
    for yr in range(low_y,hi_y):
        yl=yr*step_y;yh=yl+step_y
        domain=(u_max+yh*yh>=scale*scale)
        candidate=np.zeros(n,dtype=bool)
        for q_min,q_max,c2,i in cs:
            # h_upper <=3/5 => no transition; h_lower >=7/10 => no transition
            below=(5*scale*yh<=3*(q_min+c2*yl*yl))
            above=(10*scale*yl>=7*(q_max+c2*yh*yh))
            candidate |= ~(below|above)
        count=int(np.count_nonzero(candidate&domain))
        cover_counts.append(count)
        total_possible+=int(domain.sum())
    dx=Q(1,resolution)
    dy_contributions=[Q(n*scale*scale,2*resolution*resolution)*(Q(1,(yr*step_y)**2)-Q(1,((yr+1)*step_y)**2)) for yr,n in zip(range(low_y,hi_y),cover_counts)]
    volume=sum(dy_contributions,Q(0))*(4 if triangle else 1)
    # Domain extends below 2/3? No. Ford floor >=1/sqrt(2)>2/3.
    assert Q(2,3)**2<Q(1,2)
    core_cover=Q(73,300)
    assert volume < core_cover
    sharp_cap=Q(11,80)
    sharp_cap_proved=volume<sharp_cap
    # The published 256x128x384 run must retain its rational strict margin.
    if resolution==256 and y_den==384 and not triangle:
        assert sharp_cap_proved
    summary={"status":"EXACT_INTEGER_GEOMETRIC_SUPPORT_UPPER_BOUND","resolution":resolution,"y_den":y_den,"triangle_symmetry":triangle,
             "scale":scale,"count_of_horizontal_cells":n,"count_of_y_cells":hi_y-low_y,
             "outer_cover_voxels":sum(cover_counts),"domain_voxels":total_possible,
             "support_volume_upper_rational":f"{volume.numerator}/{volume.denominator}",
             "support_volume_upper_float_diagnostic":float(volume),
             "core_volume_upper_rational":"73/300",
             "sharp_cap_11_over_80_passed_exactly":sharp_cap_proved,
             "conditional_pointwise_residual_threshold_for_R_lt_1e-10":float(Q(1,10**10))/float(volume)**.5,
             "proof":"Exact integer interval geometry; gamma inactive when h<=3/5 or h>=7/10 throughout cuboid. If every nonidentity gamma is inactive then residual zero. All partially intersecting cubes included; hyperbolic volumes summed exactly via rationals."}
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--resolution',type=int,default=128);p.add_argument('--y-den',type=int,default=192);p.add_argument('--output',default='support_cover.json');p.add_argument('--triangle',action='store_true');args=p.parse_args()
    summary=calculate(args.resolution,args.y_den,args.triangle)
    Path(args.output).write_text(json.dumps(summary,indent=2,sort_keys=True)+'\n')
    print(json.dumps(summary,indent=2))
