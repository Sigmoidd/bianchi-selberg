"""Exact d43 finite Ford region, lower faces, face matrices and vertex isotropy."""
from fractions import Fraction as F
from itertools import product
from collections import Counter
import hashlib
import json
from pathlib import Path
from ford_exact import (add, sub, neg, mul, conj, norm, quotient, primitive,
                        candidates, center, q, find_matrix)
ROOT=Path(__file__).resolve().parent
P=((F(-1,2),F(-1,2)),(F(1,2),F(-1,2)),
   (F(1,2),F(1,2)),(F(-1,2),F(1,2)))

def affine(c,l):
    u,v=center(c,l)
    return 2*u+v,u+22*v,F(1,norm(c))-u*u-u*v-11*v*v

def clip(poly,A,B,C):
    if not poly:return ()
    out=[]
    for s,t in zip(poly,poly[1:]+poly[:1]):
        fs=A*s[0]+B*s[1]+C;ft=A*t[0]+B*t[1]+C
        if fs>=0:out.append(s)
        if (fs<0 and ft>0) or (fs>0 and ft<0):
            v=fs/(fs-ft)
            out.append((s[0]+v*(t[0]-s[0]),s[1]+v*(t[1]-s[1])))
    return tuple(dict.fromkeys(out))

def area(poly):
    return sum(s[0]*t[1]-s[1]*t[0]
               for s,t in zip(poly,poly[1:]+poly[:1]))/2

def faces():
    rows=candidates()
    if len(rows)!=1011:
        raise ArithmeticError('finite candidate enumeration changed')
    coeff=[affine(c,l) for c,l in rows]
    out=[]
    for i,(c,l) in enumerate(rows):
        ai,bi,ci=coeff[i];poly=P
        for j,(aj,bj,cj) in enumerate(coeff):
            if j==i:continue
            poly=clip(poly,ai-aj,bi-bj,ci-cj)
            if len(poly)<3:break
        if len(poly)>=3 and area(poly)>0:
            out.append((c,l,poly))
    return rows,out

def stabilizer_at(z,y2):
    """Represent PSL matrices by positive lexicographic lower-left row c."""
    cs=[(u,v) for u,v in product(range(-5,6),range(-2,3))
        if 0<norm((u,v))<=21 and (u,v)>(0,0)]
    ds=[(u,v) for u,v in product(range(-11,12),range(-4,5))
        if norm((u,v))<100]
    out=[]
    for c in cs:
        cz=mul(c,z)
        for d in ds:
            w=add(cz,d)
            if norm(w)+norm(c)*y2!=1:continue
            # z'=z on a Ford equality and determinant one iff
            # a=c*z+conj(c*z+d), with b=(a*d-1)/c.
            ra=add(cz,conj(w))
            if any(x.denominator!=1 for x in ra):continue
            a=tuple(int(x) for x in ra)
            b=quotient(sub(mul(a,d),(1,0)),c)
            if b is None:continue
            if sub(mul(a,d),mul(b,c))!=(1,0):
                raise ArithmeticError('stabilizer determinant')
            out.append((a,b,c,d))
    return out

def result():
    rows,active=faces()
    face_entries=[];vertices={}
    for c,l,poly in active:
        matrix=find_matrix(c,l)
        if sub(mul(matrix[0],matrix[3]),mul(matrix[1],matrix[2]))!=(1,0):
            raise ArithmeticError('face determinant')
        for z in poly:
            y2=q(c,l,*z)
            if y2<=0:raise ArithmeticError('nonpositive vertex height')
            old=vertices.setdefault(z,y2)
            if old!=y2:raise ArithmeticError('nonmatching vertex heights')
        face_entries.append(dict(c=c,l=l,matrix=matrix,
            polygon=[[str(a),str(b)] for a,b in poly],area=str(area(poly))))
    if sum((area(poly) for _,_,poly in active),F())!=1:
        raise ArithmeticError('floor polygons do not tile the unit cell')
    lowest=min(vertices.values())
    if lowest!=F(2,43):raise ArithmeticError('unexpected lowest vertex')
    stabilizers=[]
    for z,y2 in sorted(vertices.items()):
        mats=stabilizer_at(z,y2)
        stabilizers.append(dict(z=[str(t) for t in z],y2=str(y2),
                                order_psl=1+len(mats),matrices=mats))
    vgroups={tuple(map(F,rec['z'])):
             {tuple(tuple(x) for x in matrix) for matrix in rec['matrices']}
             for rec in stabilizers}
    edges=Counter()
    for _,_,poly in active:
        for z,w in zip(poly,poly[1:]+poly[:1]):
            edges[tuple(sorted((z,w)))]+=1
    edge_ledger=[]
    for (z,w),incidence in sorted(edges.items()):
        boundary=(z[0]==w[0] and abs(z[0])==F(1,2)) or \
                 (z[1]==w[1] and abs(z[1])==F(1,2))
        if incidence!=(1 if boundary else 2):
            raise ArithmeticError('edge incidence or vertical boundary mismatch')
        fixed=sorted(vgroups[z]&vgroups[w])
        edge_ledger.append(dict(endpoints=[[str(t) for t in p] for p in (z,w)],
                                incidence=incidence,vertical_boundary=boundary,
                                pointwise_stabilizer=fixed))
    return dict(schema='d43-finite-ford-geometry/v1',d=43,D=-43,
        subgroup='full',level_hnf=[1,0,1],
        design_sha256=hashlib.sha256((ROOT/'../../../docs/D43_TRACK_B_DESIGN.md').resolve().read_bytes()).hexdigest(),
        cover_sha256=hashlib.sha256((ROOT/'cover_witnesses.json').read_bytes()).hexdigest(),
        c_norm_max=23,lambda_norm_strict_upper=100,omitted_norm_min=25,
        cover_threshold='1/24',candidate_count=len(rows),
        active_face_count=len(active),vertex_count=len(vertices),
        edge_count=len(edge_ledger),edges=edge_ledger,
        lowest_vertex_height_squared=str(lowest),
        vertical_pairings=['T_1','T_tau'],faces=face_entries,
        vertex_stabilizers=stabilizers)

if __name__=='__main__':
    result_value=result()
    (ROOT/'geometry_result.json').write_text(json.dumps(result_value,indent=2)+'\n')
    print({k:v for k,v in result_value.items() if k not in ('faces','edges','vertex_stabilizers')})
    print('stabilizer orders',sorted({x['order_psl'] for x in result_value['vertex_stabilizers']}))
