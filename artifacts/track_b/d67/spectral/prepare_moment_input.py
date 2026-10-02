"""Generate exact integer-grid mesh input and all root-face pairing orbits."""
from fractions import Fraction as F
from itertools import permutations
from collections import defaultdict,deque
from pathlib import Path
import gzip,json,math,hashlib,argparse
from reference_mesh import ROOT
from verify_geometry import projected_map

IDENTITY=(0,1,2)
def inverse(p):return tuple(p.index(i) for i in range(3))
def compose(p,q):return tuple(p[q[i]] for i in range(3))

def prepare(directory):
    directory.mkdir(parents=True,exist_ok=True)
    raw=(ROOT/'adaptive_refinement_plan.json.gz').read_bytes();plan=json.loads(gzip.decompress(raw))
    mesh=json.loads((ROOT/'reference_mesh.json').read_text());ledger=json.loads((ROOT.parent/'geometry_result.json').read_text())
    triangles=[(tuple(map(tuple,(rec['c'],rec['l']))),tuple(tuple(map(F,z)) for z in rec['vertices'])) for rec in mesh['triangles']]
    vertices=sorted({z for _,tri in triangles for z in tri});nv=len(vertices);vi={v:i for i,v in enumerate(vertices)}
    levels=list(map(F,plan['r_levels']));layers=len(levels)-1
    zden=math.lcm(*(x.denominator for v in vertices for x in v));rden=math.lcm(*(r.denominator for r in levels))
    grid_depth=16;zs=zden*2**grid_depth;rs=rden*2**grid_depth
    if zs.bit_length()>60 or rs.bit_length()>100:raise ArithmeticError('integer-grid range')
    exact_nodes=[(*z,r) for r in levels for z in vertices]
    integers=[(int(a*zs),int(b*zs),int(r*rs)) for a,b,r in exact_nodes]
    node_id={p:i for i,p in enumerate(exact_nodes)}
    roots=[]
    for _,tri in triangles:
        a,b,c=sorted(vi[z] for z in tri)
        for k in range(layers):
            v0,v1,v2=[i+k*nv for i in (a,b,c)];w0,w1,w2=[i+(k+1)*nv for i in (a,b,c)]
            roots.extend(((v0,v1,v2,w2),(v0,v1,w1,w2),(v0,w0,w1,w2)))
    faces={tuple(sorted(t[j] for j in range(4) if j!=i)) for t in roots for i in range(4)}
    matrices={(tuple(r['c']),tuple(r['l'])):tuple(map(tuple,r['matrix'])) for r in ledger['faces']}
    floor_rows={tuple(sorted(vi[v] for v in tri)):row for row,tri in triangles}
    graph=defaultdict(list)
    def connect(src,target_ids):
        dst=tuple(sorted(target_ids));p=tuple(dst.index(n) for n in target_ids)
        if dst not in faces:raise ArithmeticError('missing paired root face')
        graph[src].append((dst,p));graph[dst].append((src,inverse(p)))
    for face in sorted(faces):
        points=[exact_nodes[n] for n in face]
        if all(p[2]==0 for p in points):
            row=floor_rows[face];matrix=matrices[row]
            mapped=[projected_map(matrix,p[:2]) for p in points]
            mean=tuple(sum((p[k] for p in mapped),F())/3 for k in (0,1));shift=tuple((x+F(1,2))//1 for x in mean)
            target=[node_id[(p[0]-shift[0],p[1]-shift[1],F())] for p in mapped]
            connect(face,target)
        for axis in (0,1):
            if all(p[axis]==F(-1,2) for p in points):
                target=[node_id[tuple(p[k]+(1 if k==axis else 0) for k in range(3))] for p in points]
                connect(face,target)
    seeds=[];visited=set();group_orders=[]
    for start in sorted(graph):
        if start in visited:continue
        frames={start:IDENTITY};queue=deque([start]);cycles=[]
        while queue:
            src=queue.popleft();ps=frames[src]
            for dst,edge in graph[src]:
                proposed=compose(edge,ps)
                if dst not in frames:frames[dst]=proposed;queue.append(dst)
                else:cycles.append(compose(inverse(frames[dst]),proposed))
        group={IDENTITY,*cycles}
        while True:
            extended={compose(p,q) for p in group for q in group}|group
            if extended==group:break
            group=extended
        if len(group)>6:raise ArithmeticError('root face permutation group')
        group_orders.append(len(group));visited.update(frames)
        for face,frame in sorted(frames.items()):seeds.append((face,start,frame,sorted(group)))
    out=directory/'moment_input.txt'
    with out.open('w') as f:
        f.write(f'{len(integers)} {len(roots)} {len(seeds)} {len(plan["leaves"])} {zs} {rs}\n')
        for p in integers:f.write(' '.join(map(str,p))+'\n')
        for t in roots:f.write(' '.join(map(str,t))+'\n')
        for face,rep,frame,group in seeds:
            f.write(' '.join(map(str,(*face,*rep,*frame,len(group),*(x for p in group for x in p))))+'\n')
        for root,path in plan['leaves']:f.write(f'{root} {path or "*"}\n')
    meta=dict(schema='d67-exact-moment-input/v1',d=67,initial_nodes=len(integers),initial_tetrahedra=len(roots),
              paired_root_faces=len(seeds),root_pairing_group_orders=sorted(set(group_orders)),leaves=len(plan['leaves']),
              z_scale=str(zs),r_scale=str(rs),grid_depth=grid_depth,
              plan_sha256=hashlib.sha256(raw).hexdigest(),input_sha256=hashlib.sha256(out.read_bytes()).hexdigest(),
              spectral_exclusion_certified=False)
    (directory/'input.json').write_text(json.dumps(meta,indent=2)+'\n')
    return meta

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('directory',type=Path);a=p.parse_args();print(json.dumps(prepare(a.directory),indent=2))
