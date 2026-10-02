"""Mapped CR threshold diagnostic. Floating quadrature: NEVER a certificate."""
from pathlib import Path
import argparse,json,sys,time
from fractions import Fraction as F
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import splu,eigsh,LinearOperator
from reference_mesh import construct,canonical,ROOT
from ford_exact import q,center
from verify_geometry import projected_map,area

G_INV=np.array([[68/67,-2/67],[-2/67,4/67]])

def planar(refinements):
    ledger=json.loads((ROOT.parent/'geometry_result.json').read_text())
    _,triangles,_,_,_=construct(ledger)
    for _ in range(refinements):
        more=[]
        for row,(a,b,c) in triangles:
            ab=tuple((a[i]+b[i])/2 for i in range(2));bc=tuple((b[i]+c[i])/2 for i in range(2));ca=tuple((c[i]+a[i])/2 for i in range(2))
            more.extend((row,canonical(t)) for t in ((a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca)))
        triangles=more
    return ledger,triangles

def build(refinements,layers):
    ledger,triangles=planar(refinements)
    vertices=sorted({v for _,tri in triangles for v in tri})
    vertex_id={v:i for i,v in enumerate(vertices)}
    # Grade g=y² exponentially at representative floor height 1/8.
    # The same r levels preserve the reference face triangulation.
    levels=[F(str(x)) for x in np.r_[0,(np.geomspace(1/8,4,layers+1)[1:]-1/8)/(4-1/8)]]
    levels[-1]=F(1)
    exact_nodes=[(a,b,r) for r in levels for a,b in vertices]
    X=np.array(exact_nodes,dtype=float)
    nv=len(vertices);tets=[];rows=[]
    for row,triangle in triangles:
        nodes=sorted(vertex_id[v] for v in triangle)
        for k in range(layers):
            v0,v1,v2=[i+k*nv for i in nodes];w0,w1,w2=[i+(k+1)*nv for i in nodes]
            tets.extend(((v0,v1,v2,w2),(v0,v1,w1,w2),(v0,w0,w1,w2)))
            rows.extend((row,)*3)
    face_id={};face_nodes=[];tet_faces=[]
    for tet in tets:
        ids=[]
        for i in range(4):
            face=tuple(sorted(tet[j] for j in range(4) if j!=i))
            if face not in face_id:
                face_id[face]=len(face_nodes);face_nodes.append(face)
            ids.append(face_id[face])
        tet_faces.append(ids)
    parent=list(range(len(face_nodes)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    def join(i,j):
        i,j=find(i),find(j)
        if i!=j:parent[max(i,j)]=min(i,j)
    node_lookup={v:i for i,v in enumerate(exact_nodes)}
    matrix_lookup={(tuple(rec['c']),tuple(rec['l'])):tuple(map(tuple,rec['matrix'])) for rec in ledger['faces']}
    bottom_rows={tuple(sorted(vertex_id[v] for v in tri)):row for row,tri in triangles}
    n_bottom=n_periodic=0
    for i,face in enumerate(face_nodes):
        nodes=[exact_nodes[j] for j in face]
        if all(v[2]==0 for v in nodes):
            row=bottom_rows[face];matrix=matrix_lookup[row]
            images=[projected_map(matrix,v[:2]) for v in nodes]
            mid=tuple(sum((v[k] for v in images),F())/3 for k in (0,1))
            shift=tuple(int((x+F(1,2))//1) for x in mid)
            target=tuple(sorted(node_lookup[(v[0]-shift[0],v[1]-shift[1],F())] for v in images))
            if target not in face_id:raise ArithmeticError('floor triangle pairing missing')
            join(i,face_id[target]);n_bottom+=1
        for axis in (0,1):
            if all(v[axis]==F(-1,2) for v in nodes):
                target=tuple(sorted(node_lookup[tuple(v[j]+(1 if j==axis else 0) for j in range(3))] for v in nodes))
                if target not in face_id:raise ArithmeticError('periodic triangle pairing missing')
                join(i,face_id[target]);n_periodic+=1
    roots=sorted({find(i) for i in range(len(parent))});lookup={r:i for i,r in enumerate(roots)}
    dofs=np.array([lookup[find(i)] for i in range(len(parent))])
    return X,tets,rows,np.array(tet_faces),face_nodes,dofs,n_bottom,n_periodic

def coefficients(row,a,b,r):
    u,v=map(float,center(*row));n=row[0][0]**2+row[0][0]*row[0][1]+17*row[0][1]**2
    da,db=a-u,b-v;H=1/n-da*da-da*db-17*db*db
    Ha=-2*da-db;Hb=-da-34*db
    D=4-H;g=(1-r)*H+4*r;h=(1-r)*np.array([Ha,Hb])/D
    mat=np.empty((3,3));mat[:2,:2]=G_INV;mat[:2,2]=-G_INV@h;mat[2,:2]=mat[:2,2];mat[2,2]=h@G_INV@h+4*g/D**2
    return D/(2*g)*mat,D/(2*g*g)

def probe(refinements=0,layers=6,rho=3):
    start=time.time();X,tets,rows,tf,face_nodes,dofs,n_bottom,n_periodic=build(refinements,layers)
    nd=int(dofs.max())+1;ii=[];jj=[];qd=[];md=[]
    qa=(5+3*np.sqrt(5))/20;qb=(5-np.sqrt(5))/20
    bary=np.full((4,4),qb);np.fill_diagonal(bary,qa)
    for tet,row,fids in zip(tets,rows,tf):
        coords=X[list(tet)];B=np.column_stack((np.ones(4),coords));inv=np.linalg.inv(B)
        volume=abs(np.linalg.det(B))/6;grads=-3*inv[1:,:]
        Q=np.zeros((4,4));M=np.zeros((4,4))
        for weights in bary:
            point=weights@coords;A,m=coefficients(row,*point);phi=1-3*weights
            Q+=volume/4*grads.T@A@grads;M+=volume/4*m*np.outer(phi,phi)
        ids=dofs[fids]
        ii.extend(np.repeat(ids,4));jj.extend(np.tile(ids,4));qd.extend(Q.ravel());md.extend(M.ravel())
    Q=coo_matrix((qd,(ii,jj)),shape=(nd,nd)).tocsc();M=coo_matrix((md,(ii,jj)),shape=(nd,nd)).tocsc()
    top=np.zeros(nd)
    for face,i in zip(face_nodes,dofs):
        if all(X[j,2]==1 for j in face):
            p=X[list(face),:2];ar=abs(np.linalg.det((p[1:]-p[0]).T))/2;top[i]+=ar
    if abs(top.sum()-1)>1e-10:raise ArithmeticError('top area mismatch')
    alpha=.5;beta=.75;K=Q+alpha*diags(top);lu=splu(K)
    mean=M@np.ones(nd);ell=mean+.25*top
    U=np.column_stack((top,ell));W=lu.solve(U)
    correction=np.linalg.inv(np.diag([-1/beta,1/rho])+U.T@W)
    def inv(x):
        y=lu.solve(x);return y-W@(correction@(U.T@y))
    def apply(x):return K@x-beta*top*(top@x)+rho*ell*(ell@x)
    op=LinearOperator((nd,nd),matvec=apply);inverse=LinearOperator((nd,nd),matvec=inv)
    values,vectors=eigsh(op,k=4,M=M,sigma=0,OPinv=inverse,tol=1e-8)
    residuals=[float(np.linalg.norm(apply(v)-mu*(M@v))/np.linalg.norm(M@v)) for mu,v in zip(values,vectors.T)]
    return dict(schema='d67-mapped-cr-threshold-diagnostic/v1',d=67,refinements=refinements,layers=layers,rho=rho,alpha=alpha,
        tetrahedra=len(tets),dofs=nd,paired_floor_faces=n_bottom,paired_translation_faces=n_periodic,
        eigenvalues=[float(x) for x in values],residuals=residuals,elapsed_seconds=time.time()-start,
        core_mass=float(mean.sum()),ell_of_constant=float(ell.sum()),penalty_functional='core_mean_plus_top_quarter',diagnostic_only=True,spectral_exclusion_certified=False,
        missing=['interval coefficient enclosures','CR error scalar bound','verified matrix positivity','independent mesh and matrix replay'])

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--refinements',type=int,default=0);p.add_argument('--layers',type=int,default=6);p.add_argument('--rho',type=float,default=3);p.add_argument('--output',type=Path);args=p.parse_args()
    result=probe(args.refinements,args.layers,args.rho)
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
