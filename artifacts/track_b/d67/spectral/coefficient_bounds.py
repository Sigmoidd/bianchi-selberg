"""Exact rational conservative pullback coefficient and CR error bounds."""
from fractions import Fraction as F
from itertools import combinations
from threshold_probe import planar
from ford_exact import q,center
from verify_geometry import area
from reference_mesh import ROOT
import argparse,json
import numpy as np

def interval_product(x,y):
    vals=[a*b for a in x for b in y];return min(vals),max(vals)

def h_bounds(row,triangle):
    u,v=center(*row)
    values=[q(*row,*z) for z in triangle];minimum=min(values)
    candidates=list(triangle)
    if all((b[0]-a[0])*(v-a[1])-(b[1]-a[1])*(u-a[0])>=0 for a,b in zip(triangle,triangle[1:]+triangle[:1])):candidates.append((u,v))
    for a,b in zip(triangle,triangle[1:]+triangle[:1]):
        dx,dy=b[0]-a[0],b[1]-a[1];da,db=a[0]-u,a[1]-v
        t=-(2*da*dx+da*dy+db*dx+34*db*dy)/(2*(dx*dx+dx*dy+17*dy*dy))
        t=max(F(),min(F(1),t));candidates.append((a[0]+t*dx,a[1]+t*dy))
    return minimum,max(q(*row,*z) for z in candidates)

def bounds(row,triangle,rlo,rhi):
    hmin,hmax=h_bounds(row,triangle);dmin,dmax=4-hmax,4-hmin
    gmin=(1-rlo)*hmin+4*rlo;gmax=(1-rhi)*hmax+4*rhi
    if hmin<=0 or hmax>=4:raise ArithmeticError('invalid reference height')
    amin=dmin/(2*gmax);bmin=2/dmax;mass=dmax/(2*gmin*gmin)
    u,v=center(*row);ha=[-2*(a-u)-(b-v) for a,b in triangle];hb=[-(a-u)-34*(b-v) for a,b in triangle]
    intervals=[]
    for vals in (ha,hb):
        p=interval_product((min(vals),max(vals)),(1-rhi,1-rlo));intervals.append(interval_product(p,(1/dmax,1/dmin)))
    mid=[(a+b)/2 for a,b in intervals];rad=[(b-a)/2 for a,b in intervals]
    R2=(68*rad[0]**2+4*rad[0]*rad[1]+4*rad[1]**2)/67
    # Reverse Young: epsilon=2 amin R²/(bmin+2 amin R²).
    alower=amin*bmin/(bmin+2*amin*R2);blower=bmin/2
    return alower,blower,tuple(mid),mass

def volume(t):
    a,b,c,d=t
    x=[tuple(v[i]-a[i] for i in range(3)) for v in (b,c,d)]
    det=x[0][0]*(x[1][1]*x[2][2]-x[1][2]*x[2][1])-x[0][1]*(x[1][0]*x[2][2]-x[1][2]*x[2][0])+x[0][2]*(x[1][0]*x[2][1]-x[1][1]*x[2][0])
    return abs(det)/6

def metric_diameter(t,alower,blower,mid):
    values=[]
    for a,b in combinations(t,2):
        dx,dy,dr=[a[i]-b[i] for i in range(3)]
        values.append((dx*dx+dx*dy+17*dy*dy)/alower+(dr+mid[0]*dx+mid[1]*dy)**2/blower)
    return max(values)

def error_bounds(refinements=0,layers=6):
    _,triangles=planar(refinements)
    levels=[F(str(x)) for x in np.r_[0,(np.geomspace(1/8,4,layers+1)[1:]-1/8)/(4-1/8)]];levels[-1]=F(1)
    gamma=F();sigma=F();worst=None;n=0
    for row,tri in triangles:
        a,b,c=sorted(tri)
        for rlo,rhi in zip(levels,levels[1:]):
            alower,blower,mid,mass=bounds(row,tri,rlo,rhi)
            v0,v1,v2=[(*p,rlo) for p in (a,b,c)];w0,w1,w2=[(*p,rhi) for p in (a,b,c)]
            for tet in ((v0,v1,v2,w2),(v0,v1,w1,w2),(v0,w0,w1,w2)):
                diameter=metric_diameter(tet,alower,blower,mid);local=mass*F(8,45)*diameter
                if local>gamma:gamma=local;worst=dict(c=row[0],l=row[1],rlo=str(rlo),rhi=str(rhi),triangle=[[str(x) for x in z] for z in tri])
                term=mass*mass*volume(tet)*F(8,45)*diameter
                scaled=term*2**40
                sigma+=F(-(-scaled.numerator//scaled.denominator),2**40);n+=1
    result=dict(schema='d67-exact-mapped-cr-error-bounds/v1',d=67,refinements=refinements,layers=layers,tetrahedra=n,
        gamma_squared=str(gamma),sigma_squared=str(sigma),gamma_squared_float=float(gamma),sigma_squared_float=float(sigma),worst=worst,
        scalar_feasibility_for_any_eta=gamma<1,spectral_exclusion_certified=False)
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--refinements',type=int,default=0);p.add_argument('--layers',type=int,default=6);args=p.parse_args()
    result=error_bounds(args.refinements,args.layers)
    (ROOT/f'error_bounds_r{args.refinements}_l{args.layers}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ('gamma_squared','sigma_squared')},indent=2))
