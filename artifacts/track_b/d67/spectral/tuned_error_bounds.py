"""Exact scalar replay with a rational vertical-energy retention fraction.

The frozen mesh and original reports are preserved. This is a coefficient /
scalar ledger producer, never a threshold-matrix positivity certificate.
"""
import argparse
from fractions import Fraction as F
import gzip
import hashlib
import json
import time
from pathlib import Path
from itertools import groupby
from adaptive_error_bounds import hull, subdivide
from coefficient_bounds import h_bounds, interval_product, metric_diameter, volume
from ford_exact import center
from reference_mesh import ROOT
from verify_adaptive_plan import require, verify_plan


def envelope(row, tet, fraction):
    poly=hull([p[:2] for p in tet])
    rlo=min(p[2] for p in tet); rhi=max(p[2] for p in tet)
    hmin,hmax=h_bounds(row,poly)
    require(0<hmin<=hmax<4,'invalid height range')
    dmin,dmax=4-hmax,4-hmin
    gmin=(1-rlo)*hmin+4*rlo;gmax=(1-rhi)*hmax+4*rhi
    amin=dmin/(2*gmax);bmin=2/dmax;mass=dmax/(2*gmin*gmin)
    u,v=center(*row)
    intervals=[]
    for values in ([-2*(a-u)-(b-v) for a,b in poly],
                   [-(a-u)-34*(b-v) for a,b in poly]):
        product=interval_product((min(values),max(values)),(1-rhi,1-rlo))
        intervals.append(interval_product(product,(1/dmax,1/dmin)))
    mid=tuple((lo+hi)/2 for lo,hi in intervals)
    ra,rb=[(hi-lo)/2 for lo,hi in intervals]
    radius_squared=(68*ra*ra+4*ra*rb+4*rb*rb)/67
    reserve=(1-fraction)*bmin
    alower=amin*reserve/(reserve+amin*radius_squared)
    blower=fraction*bmin
    return alower,blower,mid,mass


def run(fraction=F(3,4),checkpoint=None):
    require(0<fraction<1,'invalid vertical retention')
    verify_plan()
    raw=(ROOT/'adaptive_refinement_plan.json.gz').read_bytes()
    plan=json.loads(gzip.decompress(raw))
    levels=list(map(F,plan['r_levels']))
    mesh=json.loads((ROOT/'reference_mesh.json').read_text())
    initial=[]
    for rec in mesh['triangles']:
        row=tuple(map(tuple,(rec['c'],rec['l'])))
        a,b,c=sorted(tuple(map(F,p)) for p in rec['vertices'])
        for lo,hi in zip(levels,levels[1:]):
            v0,v1,v2=[(*p,lo) for p in (a,b,c)]
            w0,w1,w2=[(*p,hi) for p in (a,b,c)]
            for tet in ((v0,v1,v2,w2),(v0,v1,w1,w2),(v0,w0,w1,w2)):
                initial.append((row,tet))
    gamma=F();sigma=F();count=0;worst=None;start=time.monotonic();last=start;resume_root=0
    binding=dict(plan_sha256=hashlib.sha256(raw).hexdigest(),vertical_fraction=str(fraction),
                 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if checkpoint and checkpoint.exists():
        saved=json.loads(checkpoint.read_text())
        require(saved['binding']==binding,'scalar checkpoint binding mismatch')
        gamma=F(saved['gamma']);sigma=F(saved['sigma']);count=saved['count']
        worst=saved['worst'];resume_root=saved['next_root']
    for root,paths in groupby(plan['leaves'],key=lambda r:r[0]):
        if root<resume_root:continue
        row,base=initial[root];base_volume=volume(base)
        # Share reconstructed prefixes, never round coordinates.
        cache={'':base}
        for _,path in paths:
            for k in range(1,len(path)+1):
                prefix=path[:k]
                if prefix not in cache:
                    children=subdivide(cache[prefix[:-1]])
                    for j,t in enumerate(children):cache[prefix[:-1]+str(j)]=t
            tet=cache[path]
            a,b,mid,mass=envelope(row,tet,fraction)
            local_gamma=mass*F(1661,15000)*metric_diameter(tet,a,b,mid)
            if local_gamma>gamma:gamma=local_gamma;worst=[root,path]
            term=mass*base_volume/F(8**len(path))*local_gamma
            dyadic=term*2**40
            sigma+=F(-(-dyadic.numerator//dyadic.denominator),2**40)
            count+=1
            now=time.monotonic()
            if now-last>=30:
                print(json.dumps(dict(leaves=count,total=len(plan['leaves']),
                    gamma_squared_upper=float(gamma),elapsed_seconds=now-start)),flush=True)
                last=now
        if checkpoint:
            saved=dict(binding=binding,gamma=str(gamma),sigma=str(sigma),count=count,
                       worst=worst,next_root=root+1)
            temporary=checkpoint.with_suffix('.tmp')
            temporary.write_text(json.dumps(saved))
            temporary.replace(checkpoint)
    require(count==len(plan['leaves']),'incomplete scalar replay')
    eta=F(1,10);theta=F(9,10);rho=F(5)
    c_e=1-(1+1/eta)*gamma-rho*(1/theta-1)*sigma
    return dict(schema='d67-tuned-exact-scalar-ledger/v1',d=67,
        vertical_energy_fraction=str(fraction),kappa_squared_upper='1661/15000',
        plan_sha256=hashlib.sha256(raw).hexdigest(),leaf_tetrahedra=count,
        scalar_producer_sha256=binding['source_sha256'],
        gamma_squared_upper=str(gamma),sigma_squared_upper=str(sigma),
        gamma_squared_float=float(gamma),sigma_squared_float=float(sigma),
        eta=str(eta),theta=str(theta),rho=str(rho),c_e=str(c_e),c_e_float=float(c_e),
        worst_root_and_path=worst,scalar_budget_pass=c_e>=0,
        threshold_matrix_positivity_verified=False,spectral_exclusion_certified=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--vertical-fraction',default='3/4');p.add_argument('--output',type=Path);p.add_argument('--checkpoint',type=Path);args=p.parse_args()
    result=run(F(args.vertical_fraction),args.checkpoint);data=json.dumps(result,indent=2)+'\n'
    if args.output:args.output.write_text(data)
    print(data,end='')
