"""Exact obstruction to testing the pencil on all H1(K), not the glued core.

The trial v(z,y)=b has both mean functionals zero by central symmetry.
No sampled quadrature or approximate eigenvalue is used in this witness.
"""
from fractions import Fraction as F
from math import factorial
import json
from ford_exact import candidates,center,q
from verify_cover import verify as verify_floor_cover

def clamp(x,lo,hi):return max(lo,min(x,hi))

def maximum(row,a,b,w):
    c,l=row;u,v=center(c,l)
    points=[]
    if a<=u<=a+w and b<=v<=b+w:points.append((u,v))
    for x in (a,a+w):points.append((x,clamp(v-(x-u)/34,b,b+w)))
    for y in (b,b+w):points.append((clamp(u-(y-v)/2,a,a+w),y))
    return max(q(c,l,*z) for z in points)

def certify(n=16):
    floor=verify_floor_cover()
    rows=candidates();mass=F();w=F(1,n)
    if len(rows)!=1575:raise ArithmeticError('candidate enumeration mismatch')
    # Verify the row set is invariant under z -> -z (lambda -> -lambda).
    if {(c,(-l[0],-l[1])) for c,l in rows}!=set(rows):
        raise ArithmeticError('central symmetry missing')
    for i in range(n):
        for j in range(n):
            a=F(i,n)-F(1,2);b=F(j,n)-F(1,2)
            h=max(maximum(row,a,b,w) for row in rows)
            if not 0<h<=1:raise ArithmeticError('cell floor bound')
            mass+=w*((b+w)**3-b**3)/3*(F(1,2)/h-F(1,8))
    # q>=1/34. Q/A = (4/67) integral_P log(2/sqrt(q)) da db.
    # log(2 sqrt(34))<5/2 because exp(5)>136, proved by its series.
    exp5_lower=sum((F(5**k,factorial(k)) for k in range(12)),F())
    if not exp5_lower>136:raise ArithmeticError('log upper bound')
    energy=F(10,67)
    # s=1/2, lambda=3/4. L_s=t=0, so B_s=Q-3M/4<0.
    if not mass>F(6,25):raise ArithmeticError('simple mass bound failed')
    defect=energy-F(3,4)*F(6,25)
    if defect>=0:raise ArithmeticError('obstruction not certified; refine grid')
    return dict(schema='d67-relaxed-core-obstruction/v1',d=67,
                grid=n,candidate_rows=len(rows),trial='v=b',s='1/2',
                design_sha256=floor['design_sha256'],
                cover_sha256=floor['cover_sha256'],
                Q_over_cusp_area_strict_upper=str(energy),
                M_over_cusp_area_strict_lower='6/25',
                exact_cell_mass_lower=str(mass),
                B_over_cusp_area_strict_upper=str(defect),
                exp5_rational_lower=str(exp5_lower),
                relaxed_H1_criterion_obstructed=True,
                trial_satisfies_quotient_face_identifications=False,
                spectral_exclusion_verified=False)

if __name__=='__main__':print(json.dumps(certify(),indent=2))
