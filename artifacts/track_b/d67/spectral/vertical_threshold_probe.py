"""Floating polynomial-subspace diagnostic on the complete adaptive mesh.

This uses conservative coefficient envelopes but floating arithmetic and
quadrature. It cannot prove either exclusion or a spectral obstruction.
"""
import argparse
import json
from pathlib import Path
from itertools import combinations
import numpy as np
from scipy.linalg import eigh
from reference_mesh import ROOT
from ford_exact import center


def load(prefix):
    with prefix.with_suffix('.topology.bin').open('rb') as source:
        header = np.fromfile(source, '<u8', 6)
        if header[0] != 0x4436374d4f4d3031:
            raise ArithmeticError('topology schema')
        nn, nf, nc, nt, nd = map(int, header[1:])
        nodes = np.fromfile(source, dtype=np.dtype([
            ('a', '<i8'), ('b', '<i8'), ('lo', '<u8'), ('hi', '<i8'),
            ('parents', '<u4', (2,))]), count=nn)
        source.seek(nf * 36 + nc * 64, 1)
        leaves = np.fromfile(source, dtype=np.dtype([
            ('nodes', '<u4', (4,)), ('faces', '<u4', (4,)),
            ('root', '<u4'), ('depth', '<u4'), ('path', '<u8')]), count=nt)
        if len(leaves) != nt or source.read(1):
            raise ArithmeticError('topology layout mismatch')
    inp = json.loads((prefix.parent / 'input.json').read_text())
    x = np.column_stack((nodes['a'] / float(inp['z_scale']),
                         nodes['b'] / float(inp['z_scale']),
                         (nodes['lo'].astype(float) +
                          nodes['hi'].astype(float) * 2.0**64) / float(inp['r_scale'])))
    return x, leaves, nd


def probe(prefix, degree=6, vertical_fraction=.5):
    x, leaves, nd = load(prefix)
    mesh = json.loads((ROOT / 'reference_mesh.json').read_text())
    centers = []; inverse_norms = []
    for rec in mesh['triangles']:
        c, l = tuple(rec['c']), tuple(rec['l'])
        centers.append(list(map(float, center(c, l))))
        inverse_norms.append(1 / (c[0]**2 + c[0]*c[1] + 17*c[1]**2))
    centers = np.array(centers); inverse_norms = np.array(inverse_norms)
    dim = degree + 1
    Q = np.zeros((dim, dim)); M = Q.copy(); mean = np.zeros(dim)
    gamma=0.; sigma=0.
    qa=(5+3*np.sqrt(5))/20; qb=(5-np.sqrt(5))/20
    bary = np.full((4,4), qb); np.fill_diagonal(bary, qa)
    GI = np.array([[68,-2],[-2,4]]) / 67
    for start in range(0, len(leaves), 10000):
        batch = leaves[start:start+10000]; pts = x[batch['nodes']]
        owner = batch['root'] // 36
        z = pts[:,:,:2] - centers[owner,None,:]
        da, db = z[:,:,0], z[:,:,1]
        H = inverse_norms[owner,None] - da*da-da*db-17*db*db
        hmin = H.min(axis=1); hmax = H.max(axis=1)
        for i,j in combinations(range(4),2):
            dx,dy = da[:,j]-da[:,i], db[:,j]-db[:,i]
            den = 2*(dx*dx+dx*dy+17*dy*dy)
            numerator = -(2*da[:,i]*dx+da[:,i]*dy+db[:,i]*dx+34*db[:,i]*dy)
            t = np.clip(np.divide(numerator,den,out=np.zeros_like(den),where=den!=0),0,1)
            a,b=da[:,i]+t*dx,db[:,i]+t*dy
            hmax=np.maximum(hmax,inverse_norms[owner]-a*a-a*b-17*b*b)
        insidebox=(da.min(axis=1)<=0)&(da.max(axis=1)>=0)&(db.min(axis=1)<=0)&(db.max(axis=1)>=0)
        hmax=np.where(insidebox,inverse_norms[owner],hmax)
        dmin,dmax=4-hmax,4-hmin
        rlo,rhi=pts[:,:,2].min(axis=1),pts[:,:,2].max(axis=1)
        gmin=(1-rlo)*hmin+4*rlo; gmax=(1-rhi)*hmax+4*rhi
        amin=dmin/(2*gmax); bmin=2/dmax; mass=dmax/(2*gmin*gmin)
        bounds=[]
        for deriv in (-2*da-db,-da-34*db):
            lo,hi=deriv.min(axis=1),deriv.max(axis=1)
            products=np.stack([lo*(1-rhi),lo*(1-rlo),hi*(1-rhi),hi*(1-rlo)])
            plo,phi=products.min(axis=0),products.max(axis=0)
            vals=np.stack([plo/dmin,plo/dmax,phi/dmin,phi/dmax])
            bounds.append((vals.min(axis=0),vals.max(axis=0)))
        mid=np.column_stack([(lo+hi)/2 for lo,hi in bounds])
        rad=np.column_stack([(hi-lo)/2 for lo,hi in bounds])
        R2=(68*rad[:,0]**2+4*rad[:,0]*rad[:,1]+4*rad[:,1]**2)/67
        reserve=(1-vertical_fraction)*bmin
        alower=amin*reserve/(reserve+amin*R2); blower=vertical_fraction*bmin
        B=np.zeros((len(batch),3,3)); B[:,:2,:2]=alower[:,None,None]*GI
        off=-alower[:,None]*(mid@GI); B[:,:2,2]=off; B[:,2,:2]=off
        B[:,2,2]=alower*np.einsum('bi,ij,bj->b',mid,GI,mid)+blower
        affine=np.concatenate((np.ones((len(batch),4,1)),pts),axis=2)
        inv=np.linalg.inv(affine); volume=np.abs(np.linalg.det(affine))/6
        diameter=np.zeros(len(batch))
        for i,j in combinations(range(4),2):
            delta=pts[:,i,:]-pts[:,j,:]; dx,dy,dr=delta.T
            metric=(dx*dx+dx*dy+17*dy*dy)/alower+(dr+mid[:,0]*dx+mid[:,1]*dy)**2/blower
            diameter=np.maximum(diameter,metric)
        local_gamma=mass*(1661/15000)*diameter
        gamma=max(gamma,float(local_gamma.max()))
        sigma+=float(np.sum(mass*volume*local_gamma))
        grads=-3*inv[:,1:,:]
        # Triangle mean r^k: 2 h_k(r0,r1,r2)/((k+1)(k+2)).
        values=np.empty((len(batch),4,dim))
        for omit in range(4):
            r=pts[:,[j for j in range(4) if j!=omit],2]
            h=np.zeros((len(batch),dim)); h[:,0]=1
            for v in range(3):
                for k in range(1,dim):h[:,k]+=r[:,v]*h[:,k-1]
            values[:,omit,:]=h * np.array([2/((k+1)*(k+2)) for k in range(dim)])
        localgrad=np.einsum('bif,bfk->bik',grads,values)
        Q+=np.einsum('bik,bij,bjl,b->kl',localgrad,B,localgrad,volume)
        M+=np.einsum('bfk,bfl,b->kl',values,values,.45*volume*mass)
        total=values.sum(axis=1)
        M-=np.einsum('bk,bl,b->kl',total,total,.05*volume*mass)
        for weights in bary:
            point=np.einsum('i,bij->bj',weights,pts)
            a,b=point[:,0]-centers[owner,0],point[:,1]-centers[owner,1]
            Hq=inverse_norms[owner]-a*a-a*b-17*b*b
            mq=(4-Hq)/(2*((1-point[:,2])*Hq+4*point[:,2])**2)
            phi=1-3*weights
            uq=np.einsum('f,bfk->bk',phi,values)
            mean+=np.einsum('bk,b->k',uq,volume*mq/4)
    top=np.ones(dim); ell=mean+.25*top
    N=Q-.25*np.outer(top,top)+.5*np.outer(ell,ell)-1.1*M
    eigen, vectors=eigh(N,M)
    return dict(schema='d67-adaptive-vertical-polynomial-probe/v1',degree=degree,
                leaves=len(leaves),master_dofs=nd,minimum_restricted_eigenvalue=float(eigen[0]),
                trial_polynomial_coefficients=vectors[:,0].tolist(),
                normalized_trial_threshold_form=float(vectors[:,0]@N@vectors[:,0]),
                vertical_energy_fraction=vertical_fraction,
                gamma_squared_diagnostic=gamma,sigma_squared_diagnostic=sigma,
                scalar_c_e_diagnostic=1-11*gamma-(5/9)*sigma,
                floating_quadrature=True,coefficient_envelope='floating conservative box/hull envelope',
                spectral_exclusion_certified=False)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('prefix',type=Path);p.add_argument('--degree',type=int,default=6);p.add_argument('--vertical-fraction',type=float,default=.5);p.add_argument('--output',type=Path);args=p.parse_args()
    if not 0<args.vertical_fraction<1:raise ValueError('vertical fraction must lie in (0,1)')
    data=json.dumps(probe(args.prefix,args.degree,args.vertical_fraction),indent=2)+'\n'
    if args.output:args.output.write_text(data)
    print(data,end='')
