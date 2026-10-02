"""Independent exact replay of d67 geometry witnesses and finite face ledger."""
import json
import hashlib
from fractions import Fraction as F
from itertools import product
from collections import Counter
from functools import lru_cache
from math import gcd
from pathlib import Path
from ford_exact import (mul,add,sub,neg,conj,norm,det,center,q,primitive,quotient)
P=((F(-1,2),F(-1,2)),(F(1,2),F(-1,2)),
   (F(1,2),F(1,2)),(F(-1,2),F(1,2)))

ROOT=Path(__file__).resolve().parent

def check(ok,msg):
    if not ok: raise ArithmeticError(msg)

def pair(row):return tuple(int(x) for x in row)

def parse_polygon(v):return tuple((F(a),F(b)) for a,b in v)

def area(poly):
    return sum(x[0]*y[1]-x[1]*y[0]
               for x,y in zip(poly,poly[1:]+poly[:1]))/2

def clip(poly,A,B,C):
    if not poly:return ()
    result=[]
    for i in range(len(poly)):
        x,y=poly[i],poly[(i+1)%len(poly)]
        sx=A*x[0]+B*x[1]+C
        sy=A*y[0]+B*y[1]+C
        if sx>=0:result.append(x)
        if (sx<0<sy) or (sy<0<sx):
            result.append(((sy*x[0]-sx*y[0])/(sy-sx),
                           (sy*x[1]-sx*y[1])/(sy-sx)))
    return tuple(dict.fromkeys(result))

def vertex_group(z,y2):
    # Independent reconstruction: the fixed-point height bounds N(c)<=33,
    # while both the matrix and its inverse give N(a),N(d)<196.
    found=[]
    for u,v in product(range(-7,8),range(-1,2)):
        c=(u,v)
        if not 0<norm(c)<=33 or c<=(0,0):continue
        cz=mul(c,z)
        for h,k in product(range(-16,17),range(-3,4)):
            d=(h,k);den=norm(add(cz,d))+norm(c)*y2
            if den!=1:continue
            a0=add(cz,conj(add(cz,d)))
            if not all(x.denominator==1 for x in a0):continue
            a=pair(a0)
            b=quotient(sub(mul(a,d),(1,0)),c)
            if b is not None and sub(mul(a,d),mul(b,c))==(1,0):
                found.append((a,b,c,d))
    return found

def check_group_relations(mats):
    identity=((1,0),(0,0),(0,0),(1,0))
    members={identity,*mats}
    def canonical(g):
        c=g[2]
        if c==(0,0):
            check(g==identity or g==tuple(map(neg,identity)),
                  'nonidentity upper triangular point stabilizer')
            return identity
        return g if c>(0,0) else tuple(map(neg,g))
    for g in members:
        a,b,c,d=g
        inverse=(d,neg(b),neg(c),a)
        check(canonical(inverse) in members,'stabilizer inverse missing')
        for h in members:
            e,f,k,l=h
            prod_matrix=(add(mul(a,e),mul(b,k)),add(mul(a,f),mul(b,l)),
                         add(mul(c,e),mul(d,k)),add(mul(c,f),mul(d,l)))
            check(canonical(prod_matrix) in members,
                  'stabilizer composition missing')

def candidate_rows():
    # Reconstruct from the proved norm bounds, not from serialized faces.
    cs=[(a,b) for a,b in product(range(-6,7),range(-1,2))
        if 0<norm((a,b))<=29 and (a,b)>(0,0)]
    ls=[(a,b) for a,b in product(range(-14,15),range(-3,4))
        if norm((a,b))<169]
    return [(c,l) for c in cs for l in ls if primitive(c,l)]

@lru_cache(maxsize=None)
def affine_row(row):
    u,v=center(*row)
    return 2*u+v,u+34*v,F(1,norm(row[0]))-u*u-u*v-17*v*v

def halfspace(row,other):
    return tuple(a-b for a,b in zip(affine_row(row),affine_row(other)))

def projected_map(matrix,z):
    a,b,c,d=map(pair,matrix)
    def divide(x,y):
        p=mul(x,conj(y));n=norm(y)
        return F(p[0],n),F(p[1],n)
    return sub(divide(a,c),divide(conj(add(mul(c,z),d)),c))

def verify_cover(rows,data):
    from verify_cover import verify as replay_cover
    replay_cover(data)
    return len(data)

def verify(result=None,cover=None):
    if result is None:result=json.loads((ROOT/'geometry_result.json').read_text())
    if cover is None:cover=json.loads((ROOT/'cover_witnesses.json').read_text())
    check((result['d'],result['D'])==(67,-67),'field identity mismatch')
    check(result['schema']=='d67-finite-ford-geometry/v1','schema mismatch')
    check((result['subgroup'],result['level_hnf'])==('full',[1,0,1]),
          'full level-one identity mismatch')
    check(result['design_sha256']==hashlib.sha256(
          (ROOT/'../../../docs/D67_TRACK_B_DESIGN.md').resolve().read_bytes()).hexdigest(),
          'frozen design hash mismatch')
    check(result['cover_sha256']==hashlib.sha256(
          json.dumps(cover,indent=2).encode('utf-8')).hexdigest(),
          'cover witness hash mismatch')
    check((result['c_norm_max'],result['lambda_norm_strict_upper'],
           result['omitted_norm_min'],result['cover_threshold'])==(29,169,36,'1/34'),
          'finite termination bounds mismatch')
    check(result['vertical_pairings']==['T_1','T_tau'],'vertical pairing mismatch')
    rows=candidate_rows()
    check(len(rows)==result['candidate_count']==1575,'candidate enumeration mismatch')
    cover_n=verify_cover(rows,cover)
    active=[]
    for row in rows:
        poly=P
        for other in rows:
            if other==row:continue
            A,B,C=halfspace(row,other)
            poly=clip(poly,A,B,C)
            if len(poly)<3:break
        if len(poly)>=3 and area(poly)>0:
            active.append((row,poly))
    stored=result['faces']
    check(len(active)==len(stored)==result['active_face_count']==37,'active face mismatch')
    check(sum((area(p) for _,p in active),F())==1,'floor partition area')
    vertices={}
    for i,((c,l),poly) in enumerate(active):
        rec=stored[i]
        check((c,l)==(pair(rec['c']),pair(rec['l'])),'face row mismatch')
        check(poly==parse_polygon(rec['polygon']),'face polygon mismatch')
        check(area(poly)==F(rec['area'])>0,'face area mismatch')
        matrix=tuple(map(pair,rec['matrix']))
        a,b,cc,d=matrix
        check(cc==c and d==neg(l) and
              sub(mul(a,d),mul(b,cc))==(1,0),'face pairing determinant')
        for z in poly:
            v=q(c,l,*z)
            check(v>0,'nonpositive face vertex')
            check(vertices.setdefault(z,v)==v,'face height disagreement')
        # On its hemisphere the matrix acts as an affine real-plane isometry.
        # Clip its projected image into translation cells and check each
        # target piece against every independently enumerated inequality.
        mapped=tuple(projected_map(matrix,z) for z in poly)
        check(abs(area(mapped))==area(poly),'pairing area/isometry mismatch')
        total=F()
        for m,n in product(range(-2,3),repeat=2):
            image=clip(mapped,F(1),F(0),F(1,2)-m)
            image=clip(image,F(-1),F(0),F(1,2)+m)
            image=clip(image,F(0),F(1),F(1,2)-n)
            image=clip(image,F(0),F(-1),F(1,2)+n)
            if len(image)<3 or area(image)==0:continue
            total+=abs(area(image))
            target_l=sub(a,mul(c,(m,n)))
            for z in image:
                zz=sub(z,(m,n))
                height=q(c,target_l,*zz)
                check(height>0 and all(q(C,L,*zz)<=height for C,L in rows),
                      f'face {i} image is not on the reconstructed floor')
        check(total==area(poly),'pairing image has gaps/overlaps')
    check(len(vertices)==result['vertex_count']==66,'vertex count mismatch')
    check(min(vertices.values())==F(result['lowest_vertex_height_squared'])==F(2,67),
          'minimum height mismatch')
    stored_stab=result['vertex_stabilizers']
    check(len(stored_stab)==len(vertices),'stabilizer ledger incomplete')
    for (z,y2),rec in zip(sorted(vertices.items()),stored_stab):
        check(tuple(map(F,rec['z']))==z and F(rec['y2'])==y2,'vertex mismatch')
        mats=vertex_group(z,y2)
        check([tuple(map(pair,m)) for m in rec['matrices']]==mats,
              'vertex stabilizer mismatch')
        check(rec['order_psl']==len(mats)+1,'stabilizer order mismatch')
        check_group_relations(mats)
    groups={tuple(map(F,x['z'])):
            {tuple(map(pair,m)) for m in x['matrices']}
            for x in stored_stab}
    edge_counts=Counter()
    for _,poly in active:
        for z,w in zip(poly,poly[1:]+poly[:1]):
            edge_counts[tuple(sorted((z,w)))]+=1
    check(len(edge_counts)==result['edge_count']==len(result['edges'])==102,
          'edge ledger count mismatch')
    for ((z,w),incidence),record in zip(sorted(edge_counts.items()),result['edges']):
        check(tuple(tuple(map(F,p)) for p in record['endpoints'])==(z,w),
              'edge endpoints mismatch')
        vertical=(z[0]==w[0] and abs(z[0])==F(1,2)) or \
                 (z[1]==w[1] and abs(z[1])==F(1,2))
        check(incidence==record['incidence']==(1 if vertical else 2) and
              record['vertical_boundary']==vertical,'edge incidence mismatch')
        fixed=groups[z]&groups[w]
        check({tuple(map(pair,m)) for m in record['pointwise_stabilizer']}==fixed,
              'edge pointwise stabilizer mismatch')
    # A nonidentity orientation-preserving isometry cannot fix an open
    # two-dimensional face pointwise. Edge pointwise groups are obtained
    # as intersections of its two endpoint vertex groups.
    return dict(verified=True,candidate_rows=len(rows),cover_leaves=cover_n,
                active_faces=len(active),edges=len(edge_counts),vertices=len(vertices),
                lowest_height_squared='2/67',
                stabilizer_orders=sorted({x['order_psl'] for x in stored_stab}))

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--geometry',type=Path,default=ROOT/'geometry_result.json')
    parser.add_argument('--cover',type=Path,default=ROOT/'cover_witnesses.json')
    args=parser.parse_args()
    print(json.dumps(verify(json.loads(args.geometry.read_text()),
                            json.loads(args.cover.read_text())),indent=2))
