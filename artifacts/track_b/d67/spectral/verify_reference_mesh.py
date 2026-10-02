"""Independent exact mesh checks; no numerical certificate promotion."""
import json,hashlib,sys
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from itertools import combinations
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
from ford_exact import q,norm
from verify_geometry import candidate_rows,halfspace,projected_map,parse_polygon,area,verify,clip

def require(ok,message):
    if not ok:raise ArithmeticError(message)

def ordered(t):
    if area(t)<0:t=tuple(reversed(t))
    k=min(range(len(t)),key=lambda i:t[i]);return t[k:]+t[:k]

def overlap(t,u):
    for a,b in zip(u,u[1:]+u[:1]):
        t=clip(t,a[1]-b[1],b[0]-a[0],a[0]*b[1]-a[1]*b[0])
        if len(t)<3:return F()
    return abs(area(t))

def verify_mesh(result=None,replay_geometry=True):
    if result is None:result=json.loads((ROOT/'reference_mesh.json').read_text())
    require(result['schema']=='d67-pairing-reference-mesh/v1' and result['d']==67 and result['full_level_one'],'mesh identity')
    require(result['spectral_exclusion_certified'] is False,'unearned spectral claim')
    data=(ROOT.parent/'geometry_result.json').read_bytes()
    require(result['geometry_sha256']==hashlib.sha256(data).hexdigest(),'geometry hash')
    ledger=json.loads(data)
    if replay_geometry:verify(ledger)
    mats={(tuple(r['c']),tuple(r['l'])):tuple(map(tuple,r['matrix'])) for r in ledger['faces']}
    triangles=[];edges=Counter();rows=candidate_rows();total=F()
    for rec in result['triangles']:
        key=(tuple(rec['c']),tuple(rec['l']));t=parse_polygon(rec['vertices'])
        require(key in mats,'unknown floor row')
        require(len(t)==3 and area(t)>0,'triangle orientation')
        for z in t:
            require(all(-F(1,2)<=x<=F(1,2) for x in z),'triangle outside cell')
            require(q(*key,*z)>0,'nonpositive floor')
            for row in rows:
                A,B,C=halfspace(key,row)
                require(A*z[0]+B*z[1]+C>=0,'triangle outside active face')
        total+=area(t);triangles.append((key,t))
        for a,b in zip(t,t[1:]+t[:1]):edges[tuple(sorted((a,b)))]+=1
    require(total==1 and len(triangles)==result['triangle_count'],'area/count')
    for i,(_,t) in enumerate(triangles):
        for _,u in triangles[i+1:]:
            if max(z[0] for z in t)<=min(z[0] for z in u) or max(z[0] for z in u)<=min(z[0] for z in t) or max(z[1] for z in t)<=min(z[1] for z in u) or max(z[1] for z in u)<=min(z[1] for z in t):continue
            require(overlap(t,u)==0,'triangle interiors overlap')
    require(len({t for _,t in triangles})==len(triangles),'duplicate triangles')
    G=((F(1),F(1,2)),(F(1,2),F(17)))
    for i,(key,t) in enumerate(triangles):
        matrix=mats[key];j=result['triangles'][i]['paired_triangle']
        require(isinstance(j,int) and 0<=j<len(triangles),'pairing index')
        targetkey,target=triangles[j]
        image=tuple(projected_map(matrix,z) for z in t)
        mean=tuple(sum((z[k] for z in image),F())/3 for k in (0,1))
        shift=tuple((x+F(1,2))//1 for x in mean)
        mapped=tuple((z[0]-shift[0],z[1]-shift[1]) for z in image)
        require(ordered(mapped)==target,'entire paired triangle mismatch')
        require(abs(area(mapped))==area(t),'paired area')
        origin=projected_map(matrix,(F(),F()))
        columns=[tuple(projected_map(matrix,e)[k]-origin[k] for k in (0,1)) for e in ((F(1),F()),(F(),F(1)))]
        for a in (0,1):
            for b in (0,1):
                require(sum((columns[a][u]*G[u][v]*columns[b][v] for u in (0,1) for v in (0,1)),F())==G[a][b],'paired horizontal metric')
        # Equal heights at the vertices plus exact metric equality prove the
        # quadratic height identity over the entire triangle.
        require(all(q(*key,*z)==q(*targetkey,*zz) for z,zz in zip(t,mapped)),'paired floor height')
    boundary=[]
    for edge,count in edges.items():
        if count==2:continue
        require(count==1,'edge incidence')
        for axis in (0,1):
            if all(z[axis]==edge[0][axis] and abs(z[axis])==F(1,2) for z in edge):
                shift=1 if edge[0][axis]<0 else -1
                partner=tuple(sorted(tuple(z[k]+(shift if k==axis else 0) for k in (0,1)) for z in edge))
                require(edges[partner]==1,'periodic boundary edge');boundary.append((edge,partner));break
        else:raise ArithmeticError('unmatched interior edge')
    stored=[tuple(tuple(tuple(F(x) for x in z) for z in edge) for edge in pair) for pair in result['side_pairs']]
    require(set(stored)==set(boundary) and len(stored)==len(boundary),'periodic edge ledger')
    return dict(verified=True,d=67,triangles=len(triangles),edges=len(edges),paired_floor_triangles=len(triangles),periodic_boundary_edges=len(boundary),spectral_exclusion_certified=False)

if __name__=='__main__':
    print(json.dumps(verify_mesh(),indent=2))
