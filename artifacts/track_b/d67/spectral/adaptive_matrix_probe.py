"""Complete adaptive finite-matrix diagnostic. NEVER a spectral certificate.

Exact moment topology is converted to sparse CSR, but coefficients, assembly
and eigenvalue calculations here are floating. Rank updates remain factored.
"""
import argparse,json,subprocess,time
from pathlib import Path
from itertools import combinations
import numpy as np
from scipy.sparse import csr_matrix,diags,save_npz,load_npz
from scipy.sparse.linalg import LinearOperator,lobpcg
from reference_mesh import ROOT
from ford_exact import center
from vertical_threshold_probe import load


def coefficient_batch(pts,owners,centers,norms,fraction=.75):
    z=pts[:,:,:2]-centers[owners,None,:]
    da,db=z[:,:,0],z[:,:,1]
    H=norms[owners,None]-da*da-da*db-17*db*db
    hmin=H.min(axis=1);hmax=H.max(axis=1)
    for i,j in combinations(range(4),2):
        dx,dy=da[:,j]-da[:,i],db[:,j]-db[:,i]
        den=2*(dx*dx+dx*dy+17*dy*dy)
        num=-(2*da[:,i]*dx+da[:,i]*dy+db[:,i]*dx+34*db[:,i]*dy)
        t=np.clip(np.divide(num,den,out=np.zeros_like(den),where=den!=0),0,1)
        a,b=da[:,i]+t*dx,db[:,i]+t*dy
        hmax=np.maximum(hmax,norms[owners]-a*a-a*b-17*b*b)
    inside=(da.min(axis=1)<=0)&(da.max(axis=1)>=0)&(db.min(axis=1)<=0)&(db.max(axis=1)>=0)
    hmax=np.where(inside,norms[owners],hmax)
    dmin,dmax=4-hmax,4-hmin
    rlo,rhi=pts[:,:,2].min(axis=1),pts[:,:,2].max(axis=1)
    gmin=(1-rlo)*hmin+4*rlo;gmax=(1-rhi)*hmax+4*rhi
    amin=dmin/(2*gmax);bmin=2/dmax;mass=dmax/(2*gmin*gmin)
    intervals=[]
    for deriv in (-2*da-db,-da-34*db):
        lo,hi=deriv.min(axis=1),deriv.max(axis=1)
        products=np.stack([lo*(1-rhi),lo*(1-rlo),hi*(1-rhi),hi*(1-rlo)])
        plo,phi=products.min(axis=0),products.max(axis=0)
        vals=np.stack([plo/dmin,plo/dmax,phi/dmin,phi/dmax])
        intervals.append((vals.min(axis=0),vals.max(axis=0)))
    mid=np.column_stack([(lo+hi)/2 for lo,hi in intervals])
    rad=np.column_stack([(hi-lo)/2 for lo,hi in intervals])
    R2=(68*rad[:,0]**2+4*rad[:,0]*rad[:,1]+4*rad[:,1]**2)/67
    if fraction=='balanced':
        p=np.minimum(.5,np.ceil(np.sqrt(amin*R2/bmin)*2**20)/2**20)
        reserve=p*bmin
        aa=np.divide(amin*reserve,reserve+amin*R2,out=amin.copy(),where=(reserve+amin*R2)!=0)
        bb=(1-p)*bmin
    else:
        reserve=(1-fraction)*bmin;aa=amin*reserve/(reserve+amin*R2);bb=fraction*bmin
    GI=np.array([[68,-2],[-2,4]])/67
    B=np.zeros((len(pts),3,3));B[:,:2,:2]=aa[:,None,None]*GI
    off=-aa[:,None]*(mid@GI);B[:,:2,2]=off;B[:,2,:2]=off
    B[:,2,2]=aa*np.einsum('bi,ij,bj->b',mid,GI,mid)+bb
    affine=np.concatenate((np.ones((len(pts),4,1)),pts),axis=2)
    inv=np.linalg.inv(affine);vol=np.abs(np.linalg.det(affine))/6
    grads=-3*inv[:,1:,:]
    q=vol[:,None,None]*np.einsum('bif,bij,bjg->bfg',grads,B,grads)
    return q,mass*vol


def integration_rule(depth):
    children=np.eye(4)[None,:,:]
    for _ in range(depth):
        a,b,c,d=np.eye(4)
        ab,ac,ad,bc,bd,cd=(a+b)/2,(a+c)/2,(a+d)/2,(b+c)/2,(b+d)/2,(c+d)/2
        split=np.array([(a,ab,ac,ad),(b,ab,bc,bd),(c,ac,bc,cd),(d,ad,bd,cd),
            (ab,ac,ad,cd),(ab,ac,bc,cd),(ab,ad,bd,cd),(ab,bc,bd,cd)])
        children=np.einsum('sij,bjk->bsik',split,children).reshape(-1,4,4)
    values=1-3*children;total=values.sum(axis=1)
    grams=(np.einsum('si,sj->sij',total,total)+np.einsum('svi,svj->sij',values,values))/20
    return children,grams


def integrated_mass(pts,owners,centers,norms,rule):
    children,grams=rule
    sub=np.einsum('sij,bjk->bsik',children,pts)
    a=sub[:,:,:,0]-centers[owners,0,None,None]
    b=sub[:,:,:,1]-centers[owners,1,None,None]
    heights=norms[owners,None,None]-a*a-a*b-17*b*b
    hmin=heights.min(axis=2);rlo=sub[:,:,:,2].min(axis=2)
    mass_upper=(4-hmin)/(2*((1-rlo)*hmin+4*rlo)**2)
    affine=np.concatenate((np.ones((len(pts),4,1)),pts),axis=2)
    volume=np.abs(np.linalg.det(affine))/6
    return np.einsum('bs,sij,b->bij',mass_upper,grams,volume/len(children))


def prolongation(prefix):
    exe=prefix.parent/'export_moment_csr'
    subprocess.run(['g++','-O2','-std=c++17',str(ROOT/'export_moment_csr.cpp'),'-o',str(exe)],check=True)
    output=prefix.with_suffix('.csr.bin')
    subprocess.run([str(exe),str(prefix.with_suffix('.rows.bin')),str(output)],check=True)
    with output.open('rb') as source:
        nr,nd,nz=map(int,np.fromfile(source,'<u8',3))
        ptr=np.fromfile(source,'<u8',nr+1).astype(np.int64)
        col=np.fromfile(source,'<u4',nz).astype(np.int32)
        val=np.fromfile(source,'<f8',nz)
    return csr_matrix((val,col,ptr),shape=(nr,nd))


def master_points(prefix,x):
    with prefix.with_suffix('.topology.bin').open('rb') as f:
        _,nn,nf,nc,_,nd=map(int,np.fromfile(f,'<u8',6))
    faces=np.memmap(prefix.with_suffix('.topology.bin'),mode='r',offset=48+nn*40,
        dtype=np.dtype([('v','<u4',(3,)),('parent','<u4'),('children','<u4',(4,)),('cell','<u4')]),shape=(nf,))
    cells=np.memmap(prefix.with_suffix('.topology.bin'),mode='r',offset=48+nn*40+nf*36,
        dtype=np.dtype([('root','<u4'),('depth','<u4'),('bary','<u4',(3,3)),('children','<u4',(4,)),('dof','<u4')]),shape=(nc,))
    selected=cells[cells['dof']!=2**32-1]
    coords=np.empty((nd,3))
    for start in range(0,len(selected),100000):
        part=selected[start:start+100000]
        weights=part['bary'].sum(axis=1).astype(float)/(3*2.0**part['depth'][:,None])
        coords[part['dof']]=np.einsum('bi,bij->bj',weights,x[faces['v'][part['root']]])
    return coords


def assemble(prefix,balanced=False,mass_depth=0):
    start=time.monotonic();x,leaves,nd=load(prefix);nt=len(leaves)
    mesh=json.loads((ROOT/'reference_mesh.json').read_text())
    centers=np.array([list(map(float,center(tuple(r['c']),tuple(r['l'])))) for r in mesh['triangles']])
    norms=np.array([1/(r['c'][0]**2+r['c'][0]*r['c'][1]+17*r['c'][1]**2) for r in mesh['triangles']])
    P=prolongation(prefix)
    qdata=np.empty((nt,4,4));mass=np.empty(nt);top=np.zeros((nt,4))
    mdata=np.empty((nt,4,4)) if mass_depth else None
    rule=integration_rule(mass_depth)
    for i in range(0,nt,10000):
        batch=leaves[i:i+10000];pts=x[batch['nodes']]
        qdata[i:i+len(batch)],mass[i:i+len(batch)]=coefficient_batch(pts,batch['root']//36,centers,norms,'balanced' if balanced else .75)
        if mass_depth:mdata[i:i+len(batch)]=integrated_mass(pts,batch['root']//36,centers,norms,rule)
        for omit in range(4):
            tri=pts[:,[k for k in range(4) if k!=omit],:]
            is_top=np.all(tri[:,:,2]==1,axis=1)
            da=tri[:,1,:2]-tri[:,0,:2];db=tri[:,2,:2]-tri[:,0,:2]
            area=np.abs(da[:,0]*db[:,1]-da[:,1]*db[:,0])/2
            top[i:i+len(batch),omit]=is_top*area
        if i%200000==0:print(json.dumps(dict(assembled_local=i,elapsed_seconds=time.monotonic()-start)),flush=True)
    indptr=np.arange(0,16*nt+1,4,dtype=np.int64)
    indices=np.broadcast_to(4*np.arange(nt,dtype=np.int32)[:,None,None]+np.arange(4,dtype=np.int32)[None,None,:],(nt,4,4)).reshape(-1).copy()
    Qlocal=csr_matrix((qdata.reshape(-1),indices,indptr),shape=(4*nt,4*nt))
    print('assembling global Q',flush=True)
    Q=(P.T@(Qlocal@P)).tocsr();del Qlocal,qdata
    if mdata is None:mdata=mass[:,None,None]*(.45*np.eye(4)-.05*np.ones((4,4)))
    Mlocal=csr_matrix((mdata.reshape(-1),indices,indptr),shape=(4*nt,4*nt))
    print('assembling global M',flush=True)
    M=(P.T@(Mlocal@P)).tocsr();del Mlocal,mdata,indices,indptr
    t=np.asarray(P.T@top.reshape(-1)).ravel()
    B=(P.T@(diags(top.reshape(-1))@P)).tocsr()
    K=(Q+.5*B).tocsr();del Q,B,P
    K.eliminate_zeros();M.eliminate_zeros()
    z=1.1*np.asarray(M@np.ones(nd)).ravel()+.25*t
    coords=master_points(prefix,x)
    np.savez(prefix.parent/'matrix_vectors.npz',t=t,z=z,coords=coords)
    save_npz(prefix.parent/'matrix_K.npz',K);save_npz(prefix.parent/'matrix_M.npz',M)
    report=dict(schema='d67-full-adaptive-floating-matrix/v1',d=67,master_dofs=nd,
        leaves=nt,coefficient_envelope='balanced' if balanced else 'fixed_three_quarters',
        finite_mass_integration_depth=mass_depth,
        mass_upper_volume=float(np.ones(nd)@(M@np.ones(nd))),
        stabilization='z=11*M_h*1/10+t_h/4; alpha=1/2',
        energy_nonzeros=K.nnz,mass_nonzeros=M.nnz,top_area=float(t.sum()),
        symmetry_defect_K=float(abs(K-K.T).max()),symmetry_defect_M=float(abs(M-M.T).max()),
        constant_energy=float(np.ones(nd)@(K@np.ones(nd))),
        floating_arithmetic=True,matrix_positivity_verified=False,spectral_exclusion_certified=False)
    (prefix.parent/'matrix_assembly.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)


def eigenprobe(directory,iterations=100,initial_vectors=None):
    import pyamg
    K=load_npz(directory/'matrix_K.npz');M=load_npz(directory/'matrix_M.npz')
    vectors=np.load(directory/'matrix_vectors.npz');t,z,x=vectors['t'],vectors['z'],vectors['coords'];nd=K.shape[0]
    def apply(v):
        if v.ndim==1:return K@v-1.1*(M@v)-.75*t*(t@v)+.5*z*(z@v)
        return K@v-1.1*(M@v)-.75*t[:,None]*(t@v)[None,:]+.5*z[:,None]*(z@v)[None,:]
    A=LinearOperator((nd,nd),matvec=apply,matmat=apply,dtype=float)
    print('building algebraic multigrid preconditioner',flush=True)
    mg=pyamg.smoothed_aggregation_solver(K,max_coarse=500,symmetry='symmetric')
    print(str(mg),flush=True)
    initial=np.column_stack((np.ones(nd),x[:,2],x[:,2]**2,
        np.sin(2*np.pi*x[:,0]),np.cos(2*np.pi*x[:,0]),
        np.sin(2*np.pi*x[:,1]),np.cos(2*np.pi*x[:,1])))
    if initial_vectors:initial=np.load(initial_vectors)
    vals,vecs,history,residuals=lobpcg(A,initial,B=M,M=mg.aspreconditioner(),largest=False,
        tol=1e-6,maxiter=iterations,retLambdaHistory=True,retResidualNormsHistory=True,verbosityLevel=1)
    np.save(directory/'diagnostic_eigenvectors.npy',vecs)
    report=dict(schema='d67-full-adaptive-floating-eigenprobe/v1',d=67,master_dofs=nd,
        eigenvalues=list(map(float,vals)),residuals=list(map(float,residuals[-1])),
        iterations=len(history),floating_arithmetic=True,matrix_positivity_verified=False,
        spectral_exclusion_certified=False)
    (directory/'matrix_eigenprobe.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('prefix',type=Path);p.add_argument('--eigenprobe',action='store_true');p.add_argument('--balanced',action='store_true');p.add_argument('--mass-depth',type=int,default=0);p.add_argument('--iterations',type=int,default=100);p.add_argument('--initial-vectors',type=Path);args=p.parse_args()
    if not 0<=args.mass_depth<=3:raise ValueError('diagnostic mass depth must lie between 0 and 3')
    if args.eigenprobe:eigenprobe(args.prefix.parent,args.iterations,args.initial_vectors)
    else:assemble(args.prefix,args.balanced,args.mass_depth)
