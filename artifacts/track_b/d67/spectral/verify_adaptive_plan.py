"""Exact adaptive refinement coverage check; not a hanging-face verifier."""
import gzip,hashlib,json
from fractions import Fraction as F
from itertools import groupby
from reference_mesh import ROOT
from verify_reference_mesh import verify_mesh

def require(ok,message):
    if not ok:raise ArithmeticError(message)

def verify_plan():
    verify_mesh()
    path=ROOT/'adaptive_refinement_plan.json.gz';raw=path.read_bytes();plan=json.loads(gzip.decompress(raw))
    report=json.loads((ROOT/'adaptive_error_bounds.json').read_text())
    require(plan['schema']=='d67-adaptive-refinement-plan/v1' and plan['d']==67,'plan identity')
    require(plan['spectral_exclusion_certified'] is False and plan['global_face_constraints_verified'] is False,'unearned plan claim')
    require(report['refinement_plan_gzip_sha256']==hashlib.sha256(raw).hexdigest(),'plan hash')
    require(plan['geometry_sha256']==hashlib.sha256((ROOT.parent/'geometry_result.json').read_bytes()).hexdigest(),'plan geometry binding')
    require(plan['reference_mesh_sha256']==hashlib.sha256((ROOT/'reference_mesh.json').read_bytes()).hexdigest(),'plan reference binding')
    levels=list(map(F,plan['r_levels']))
    require(levels[0]==0 and levels[-1]==1 and all(a<b for a,b in zip(levels,levels[1:])),'r layer partition')
    require(len(levels)==plan['layers']+1,'r layer count')
    leaves=plan['leaves'];require(leaves==sorted(leaves),'leaf ordering')
    roots=228*plan['layers']*3;seen=0;deepest=0
    for root,group in groupby(leaves,key=lambda x:x[0]):
        require(isinstance(root,int) and not isinstance(root,bool) and root==seen,'missing or invalid root tetrahedron');seen+=1;total=F();previous=None
        for _,path in group:
            require(isinstance(path,str) and all(c in '01234567' for c in path),'invalid subdivision path')
            require(previous is None or not path.startswith(previous),'overlapping/duplicate subdivision paths')
            previous=path;deepest=max(deepest,len(path));total+=F(1,8**len(path))
        require(total==1,'refined root has holes')
    require(seen==roots==report['initial_tetrahedra'],'initial root count')
    require(len(leaves)==report['leaf_tetrahedra'],'adaptive leaf count')
    require(deepest<=report['max_depth'],'adaptive depth limit')
    # Red children are an affine image of four disjoint corner tetrahedra
    # and an octahedron split into four tetrahedra around one diagonal.
    # Each has exactly 1/8 parent volume. Prefix freedom plus unit Kraft
    # sum proves exhaustive, nonoverlapping recursive coverage per root.
    return dict(verified=True,d=67,initial_tetrahedra=roots,leaf_tetrahedra=len(leaves),maximum_depth=deepest,
                coverage_verified=True,hanging_face_moment_map_verified=False,spectral_exclusion_certified=False)

if __name__=='__main__':print(json.dumps(verify_plan(),indent=2))
