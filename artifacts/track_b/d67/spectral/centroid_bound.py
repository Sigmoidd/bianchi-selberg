"""Exact moment identities and sharpened scalar budget; no spectral claim."""
from fractions import Fraction as F
from itertools import combinations
from reference_mesh import ROOT
import json,hashlib

KAPPA2=F(1661,15000)
RATIO=KAPPA2/F(8,45)

def moment_identity(vertices):
    center=tuple(sum((p[k] for p in vertices),F())/4 for k in range(3))
    # E(lambda_i lambda_j)=(1+delta_ij)/20 on a tetrahedron.
    second=sum(((F(1,10) if i==j else F(1,20))*sum((vertices[i][k]-center[k])*(vertices[j][k]-center[k]) for k in range(3)) for i in range(4) for j in range(4)),F())
    edges=sum((sum((p[k]-q[k])**2 for k in range(3)) for p,q in combinations(vertices,2)),F())
    if second!=edges/80:raise ArithmeticError('centroid second moment identity')
    return second,edges

def verify():
    witnesses=[((F(),F(),F()),(F(1),F(),F()),(F(),F(1),F()),(F(),F(),F(1))),
               ((F(-2),F(3),F(1)),(F(5),F(-1),F(2)),(F(1),F(7),F(-3)),(F(2),F(1),F(9)))]
    for t in witnesses:moment_identity(t)
    if RATIO!=F(4983,8000):raise ArithmeticError('interpolation ratio')
    base=(ROOT/'adaptive_error_bounds.json').read_bytes();r=json.loads(base)
    if r['d']!=67 or r['spectral_exclusion_certified'] is not False:raise ArithmeticError('bound scope')
    gamma=RATIO*F(r['gamma_squared_upper']);sigma=RATIO*F(r['sigma_squared_upper'])
    eta=F(1,10);theta=F(9,10);rho=F(5)
    ce=1-(1+1/eta)*gamma-rho*(1/theta-1)*sigma
    return dict(schema='d67-centroid-cr-bound/v1',d=67,kappa_squared_upper=str(KAPPA2),ratio=str(RATIO),
        source_bound_sha256=hashlib.sha256(base).hexdigest(),gamma_squared_upper=str(gamma),sigma_squared_upper=str(sigma),
        gamma_squared_float=float(gamma),sigma_squared_float=float(sigma),eta=str(eta),theta=str(theta),rho=str(rho),
        c_e=str(ce),c_e_float=float(ce),local_scalar_budget_arithmetic_pass=ce>=0,spectral_exclusion_certified=False)

if __name__=='__main__':
    print(json.dumps(verify(),indent=2))
