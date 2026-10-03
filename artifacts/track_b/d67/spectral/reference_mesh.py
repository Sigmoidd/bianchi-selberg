"""Exact reference mesh with affine Ford-floor pairings; no spectral claim."""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import sys,json,hashlib
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent))
from ford_exact import q,sub,mul,neg,norm
from verify_geometry import clip,area,parse_polygon,projected_map,verify


def canonical(poly):
    if area(poly)<0:poly=tuple(reversed(poly))
    k=min(range(len(poly)),key=lambda i:poly[i])
    return poly[k:]+poly[:k]


def intersection(poly,other):
    other=canonical(other)
    for (a,b),(c,d) in zip(other,other[1:]+other[:1]):
        poly=clip(poly,b-d,c-a,a*d-b*c)
        if len(poly)<3:return ()
    return canonical(poly) if area(poly)>0 else ()


def invert(matrix):
    a,b,c,d=matrix
    return d,neg(b),neg(c),a


def keyrow(c,l):
    return (c,l) if c>(0,0) else (neg(c),neg(l))


def target_pieces(poly,matrix):
    mapped=tuple(projected_map(matrix,z) for z in poly)
    for m,n in product(range(-2,3),repeat=2):
        image=clip(mapped,F(1),F(0),F(1,2)-m)
        image=clip(image,F(-1),F(0),F(1,2)+m)
        image=clip(image,F(0),F(1),F(1,2)-n)
        image=clip(image,F(0),F(-1),F(1,2)+n)
        if len(image)>=3 and abs(area(image))>0:
            yield (m,n),canonical(tuple(sub(z,(m,n)) for z in image))


def enrich(poly,vertices):
    result=[]
    for p,z in zip(poly,poly[1:]+poly[:1]):
        dx,dy=z[0]-p[0],z[1]-p[1]
        points=[]
        for v in vertices:
            if dx*(v[1]-p[1])!=dy*(v[0]-p[0]):continue
            t=(v[0]-p[0])/dx if dx else (v[1]-p[1])/dy
            if 0<=t<1:points.append((t,v))
        result.extend(v for t,v in sorted(points))
    return canonical(tuple(result))


def construct(ledger,rounds=12):
    records=ledger['faces']
    mats={keyrow(tuple(rec['c']),tuple(rec['l'])):tuple(map(tuple,rec['matrix'])) for rec in records}
    pieces={key:[] for key in mats}
    for rec in records:
        pieces[keyrow(tuple(rec['c']),tuple(rec['l']))]=[canonical(parse_polygon(rec['polygon']))]
    for iteration in range(rounds):
        refined={key:set() for key in mats}
        transferred_vertices=set()
        for key,polys in pieces.items():
            matrix=mats[key];a,b,c,d=matrix
            for poly in polys:
                total=F()
                for shift,image in target_pieces(poly,matrix):
                    target=keyrow(c,sub(a,mul(c,shift)))
                    if target not in pieces:raise ArithmeticError('missing image face')
                    for other in pieces[target]:
                        transferred_vertices.update(projected_map(invert(matrix),(z[0]+shift[0],z[1]+shift[1])) for z in other)
                        overlap=intersection(image,other)
                        if not overlap:continue
                        back=canonical(tuple(projected_map(invert(matrix),(z[0]+shift[0],z[1]+shift[1])) for z in overlap))
                        refined[key].add(back);total+=area(back)
                if total!=area(poly):raise ArithmeticError('pairing refinement lost area')
        all_vertices={v for polys in refined.values() for p in polys for v in p}
        all_vertices.update(transferred_vertices)
        all_vertices.update((z[0]+m,z[1]+n) for z in list(all_vertices) for m,n in product(range(-1,2),repeat=2) if -F(1,2)<=z[0]+m<=F(1,2) and -F(1,2)<=z[1]+n<=F(1,2))
        updated={key:sorted({enrich(p,all_vertices) for p in polys}) for key,polys in refined.items()}
        if updated==pieces:
            break
        pieces=updated
    else:raise ArithmeticError('pairing refinement did not terminate within cap')
    # Triangulate by affine centroids; all paired cells must then map bijectively.
    triangles=[]
    for key,polys in pieces.items():
        for poly in polys:
            centroid=tuple(sum((z[k] for z in poly),F())/len(poly) for k in range(2))
            for a,b in zip(poly,poly[1:]+poly[:1]):
                triangle=canonical((a,b,centroid))
                if area(triangle)<=0:raise ArithmeticError('degenerate triangle')
                triangles.append((key,triangle))
    lookup={t:i for i,(key,t) in enumerate(triangles)}
    pairings=[]
    for i,(key,triangle) in enumerate(triangles):
        images=list(target_pieces(triangle,mats[key]))
        if len(images)!=1:raise ArithmeticError('triangle crosses pairing translation cell')
        shift,image=images[0]
        if image not in lookup:raise ArithmeticError(f'paired floor triangulation mismatch: {i}')
        j=lookup[image]
        targetkey=triangles[j][0]
        for z,zz in zip(triangle,tuple(projected_map(mats[key],x) for x in triangle)):
            tt=sub(zz,shift)
            if q(*key,*z)!=q(*targetkey,*tt):raise ArithmeticError('paired heights mismatch')
        pairings.append(j)
    # Periodic edge agreement in both lattice directions.
    edge_counts={}
    for _,tri in triangles:
        for p,z in zip(tri,tri[1:]+tri[:1]):
            edge=tuple(sorted((p,z)));edge_counts[edge]=edge_counts.get(edge,0)+1
    side_pairs=[]
    for edge,count in edge_counts.items():
        if count==2:continue
        if count!=1:raise ArithmeticError('reference edge incidence')
        for axis in (0,1):
            if all(z[axis] in (F(-1,2),F(1,2)) for z in edge) and edge[0][axis]==edge[1][axis]:
                shift=[0,0];shift[axis]=-1 if edge[0][axis]>0 else 1
                target=tuple(sorted((z[0]+shift[0],z[1]+shift[1]) for z in edge))
                if edge_counts.get(target)!=1:raise ArithmeticError('periodic edge mesh mismatch')
                side_pairs.append((edge,target));break
        else:raise ArithmeticError('unmatched interior edge')
    if sum((area(t) for _,t in triangles),F())!=1:raise ArithmeticError('reference mesh area')
    return pieces,triangles,pairings,side_pairs,iteration+1


def result():
    ledger=json.loads((ROOT.parent/'geometry_result.json').read_text())
    pieces,triangles,pairs,sides,rounds=construct(ledger)
    return dict(schema='d67-pairing-reference-mesh/v1',d=67,full_level_one=True,
        geometry_sha256=hashlib.sha256((ROOT.parent/'geometry_result.json').read_bytes()).hexdigest(),
        refinement_rounds=rounds,cell_count=sum(map(len,pieces.values())),triangle_count=len(triangles),
        triangles=[dict(c=key[0],l=key[1],vertices=[[str(x) for x in z] for z in tri],paired_triangle=pairs[i]) for i,(key,tri) in enumerate(triangles)],
        side_pairs=[[[[str(x) for x in z] for z in edge] for edge in p] for p in sides],
        spectral_exclusion_certified=False)

if __name__=='__main__':
    r=result();(ROOT/'reference_mesh.json').write_text(json.dumps(r,indent=2)+'\n')
    print({k:v for k,v in r.items() if k not in ('triangles','side_pairs')})
