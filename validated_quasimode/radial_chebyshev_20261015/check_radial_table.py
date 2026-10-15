#!/usr/bin/env python3
"""Exact-decimal width audit plus nonrigorous midpoint cross-check of interval table."""
from decimal import Decimal,getcontext
from collections import defaultdict
from pathlib import Path
import numpy as np,mpmath as mp
from numpy.polynomial.chebyshev import chebval
getcontext().prec=120
BASE=Path(__file__).parent
rows=defaultdict(list)
for text in (BASE/'radial_coeffs_all.csv').read_text().splitlines():
 if text.startswith('#'):continue
 n,k,fl,fh,dl,dh=text.split(',')
 rows[int(n)].append((int(k),Decimal(fl),Decimal(fh),Decimal(dl),Decimal(dh)))
assert len(rows)==24 and all(len(v)==121 for v in rows.values())
worst_wf=Decimal(0);worst_wd=Decimal(0)
for n,vals in rows.items():
 assert [r[0] for r in vals]==list(range(121))
 widths_f=sum((r[2]-r[1] for r in vals),Decimal(0))
 widths_d=sum((r[4]-r[3] for r in vals),Decimal(0))
 assert all(r[1]<=r[2] and r[3]<=r[4] for r in vals)
 assert widths_f <Decimal('1e-20') and widths_d<Decimal('2e-20')
 worst_wf=max(worst_wf,widths_f);worst_wd=max(worst_wd,widths_d)
print('EXACT DECIMAL AUDIT: 24 norms × 121 coefficient intervals, nonempty and ordered')
print('MAX coefficient-width sum f =',worst_wf)
print('MAX coefficient-width sum fprime =',worst_wd)
mp.mp.dps=35
r=mp.mpf('6.622119348586011')
max_err_f=max_err_d=0.
for n in (1,5,25,50):
 v=rows[n];cf=np.array([float((x[1]+x[2])/2) for x in v]);cd=np.array([float((x[3]+x[4])/2) for x in v]);
 for y in (.6,.733,1.003,1.37,5/3):
  t=(15./8.)*(y-17./15.);z=2*mp.pi*mp.sqrt(n)*mp.mpf(str(y))
  kval=mp.re(mp.besselk(1j*r,z));kprime=-mp.re((mp.besselk(1j*r-1,z)+mp.besselk(1j*r+1,z))/2)
  target_f=float(mp.mpf(str(y))*kval);target_d=float(kval+z*kprime)
  max_err_f=max(max_err_f,abs(chebval(t,cf)-target_f))
  max_err_d=max(max_err_d,abs(chebval(t,cd)-target_d))
print('NUMERICAL ONLY: max held-out midpoint f discrepancy',max_err_f)
print('NUMERICAL ONLY: max held-out midpoint fprime discrepancy',max_err_d)
print('STATUS: nodal intervals plus mathematical ellipse/Gauss audits established; 3D residual supremum still open')
