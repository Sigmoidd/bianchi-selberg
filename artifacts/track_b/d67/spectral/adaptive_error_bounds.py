"""Adaptive exact error-bound diagnostic; no matrix or spectral certification.

Each accepted tetrahedron lies inside one original patch prism. Refinement
is local and may create hanging faces. No claim of a validated global CR
space is made: nested face-moment constraints still require reconstruction.
"""
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
import argparse,json,time,heapq,gzip,hashlib
import numpy as np
from coefficient_bounds import bounds,volume,metric_diameter
from threshold_probe import build
from reference_mesh import ROOT


def hull(points):
    points=sorted(set(points))
    def cross(o,a,b):return (a[0]-o[0])*(b[1]-o[1])-(a[1]-o[1])*(b[0]-o[0])
    lower=[]
    for p in points:
        while len(lower)>=2 and cross(lower[-2],lower[-1],p)<=0:lower.pop()
        lower.append(p)
    upper=[]
    for p in reversed(points):
        while len(upper)>=2 and cross(upper[-2],upper[-1],p)<=0:upper.pop()
        upper.append(p)
    return tuple(lower[:-1]+upper[:-1])


def subdivide(tet):
    a,b,c,d=tet
    def mid(p,q):return tuple((x+y)/2 for x,y in zip(p,q))
    ab,ac,ad,bc,bd,cd=[mid(*p) for p in ((a,b),(a,c),(a,d),(b,c),(b,d),(c,d))]
    return ((a,ab,ac,ad),(b,ab,bc,bd),(c,ac,bc,cd),(d,ad,bd,cd),
            (ab,ac,ad,cd),(ab,ac,bc,cd),(ab,ad,bd,cd),(ab,bc,bd,cd))


def local(row,tet):
    poly=hull([p[:2] for p in tet]);rlo=min(p[2] for p in tet);rhi=max(p[2] for p in tet)
    alower,blower,mid,mass=bounds(row,poly,rlo,rhi)
    diameter=metric_diameter(tet,alower,blower,mid)
    gamma=mass*F(8,45)*diameter
    return gamma,mass,volume(tet),alower,blower,mid


def run(target=F(1,20),layers=12,max_leaves=100000,max_depth=8):
    start=time.time()
    # Same rational r levels as the prototype, but generated directly.
    from threshold_probe import planar
    _,triangles=planar(0)
    levels=[F(str(x)) for x in np.r_[0,(np.geomspace(1/8,4,layers+1)[1:]-1/8)/(4-1/8)]];levels[-1]=F(1)
    queue=[]
    for row,tri in triangles:
        a,b,c=sorted(tri)
        for rlo,rhi in zip(levels,levels[1:]):
            v0,v1,v2=[(*p,rlo) for p in (a,b,c)];w0,w1,w2=[(*p,rhi) for p in (a,b,c)]
            queue.extend((row,t,0) for t in ((v0,v1,v2,w2),(v0,v1,w1,w2),(v0,w0,w1,w2)))
    initial=len(queue);heap=[];sequence=0;split=0
    for initial_id,(row,tet,depth) in enumerate(queue):
        data=local(row,tet)
        heapq.heappush(heap,(-data[0],sequence,row,tet,depth,initial_id,'',data));sequence+=1
    del queue
    while heap and -heap[0][0]>target and len(heap)+7<=max_leaves:
        _,_,row,tet,depth,initial_id,path,data=heapq.heappop(heap)
        if depth>=max_depth:raise ArithmeticError('adaptive maximum depth reached')
        for child_id,t in enumerate(subdivide(tet)):
            child_data=local(row,t)
            heapq.heappush(heap,(-child_data[0],sequence,row,t,depth+1,initial_id,path+str(child_id),child_data));sequence+=1
        split+=1
        if split%2000==0:print(json.dumps(dict(leaves=len(heap),gamma_squared_upper=float(-heap[0][0]),elapsed=time.time()-start)),flush=True)
    accepted=len(heap);depths={};gamma=F();sigma=F();unresolved=0;failures=[];plan=[]
    for _,_,row,tet,depth,initial_id,path,data in heap:
        g,mass,vol,alower,blower,mid=data
        if g>target:
            unresolved+=1
            if len(failures)<5:failures.append(dict(c=row[0],l=row[1],depth=depth,gamma_squared_upper=float(g)))
        depths[depth]=depths.get(depth,0)+1;gamma=max(gamma,g)
        term=mass*vol*g;scaled=term*2**40
        sigma+=F(-(-scaled.numerator//scaled.denominator),2**40)
        plan.append([initial_id,path])
    plan_record=dict(schema='d67-adaptive-refinement-plan/v1',d=67,layers=layers,r_levels=list(map(str,levels)),leaves=sorted(plan),global_face_constraints_verified=False,spectral_exclusion_certified=False,geometry_sha256=hashlib.sha256((ROOT.parent/'geometry_result.json').read_bytes()).hexdigest(),reference_mesh_sha256=hashlib.sha256((ROOT/'reference_mesh.json').read_bytes()).hexdigest())
    plan_bytes=gzip.compress(json.dumps(plan_record,separators=(',',':')).encode(),mtime=0)
    result=dict(schema='d67-adaptive-error-bound-diagnostic/v1',d=67,target_gamma_squared=str(target),layers=layers,
        max_leaves=max_leaves,max_depth=max_depth,initial_tetrahedra=initial,leaf_tetrahedra=accepted,subdivisions=split,depth_counts=depths,
        gamma_squared_upper=str(gamma),sigma_squared_upper=str(sigma),gamma_squared_float=float(gamma),sigma_squared_float=float(sigma),
        unresolved=unresolved,example_failures=failures,elapsed_seconds=time.time()-start,
        adaptive_face_moment_mesh_verified=False,matrix_positivity_verified=False,spectral_exclusion_certified=False)
    result['refinement_plan_gzip_sha256']=hashlib.sha256(plan_bytes).hexdigest()
    return result,plan_bytes

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--target',default='1/20');p.add_argument('--layers',type=int,default=12);p.add_argument('--max-leaves',type=int,default=100000);p.add_argument('--output',type=Path);args=p.parse_args()
    r,plan_bytes=run(F(args.target),args.layers,args.max_leaves)
    if args.output:
        args.output.write_text(json.dumps(r,indent=2)+'\n')
        args.output.with_suffix('.refinement.json.gz').write_bytes(plan_bytes)
    print(json.dumps({k:v for k,v in r.items() if k not in ('gamma_squared_upper','sigma_squared_upper')},indent=2))
