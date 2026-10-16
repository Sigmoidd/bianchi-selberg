#!/usr/bin/env python3
"""Exact-rational proof-gate audit and independent frozen 21-term reconstruction.
No floating point, no sampling, and no MPFR. It proves the arithmetic
implications only, given the explicitly stated analytic estimates.
"""
from fractions import Fraction as Q
from pathlib import Path
import json,csv,hashlib
from decimal import Decimal, getcontext
from math import isqrt
R=Path(__file__).parent
assert Q(6622119348586011,10**15) == Q('6.622119348586011')
raw=(R/'analytic_quasimode_coefficients_50_reflection_even.json').read_bytes()
d=json.loads(raw)
assert d['r']=='6.622119348586011' and d['M']==50
assert len(d['orbits'])==40
mp={}
for orb in d['orbits']:
 a=Q(orb['coefficient_real']);assert Q(orb['coefficient_imag'])==0
 assert len(orb['modes'])==4
 for sign,ab in zip((1,-1,1,-1),orb['modes']):
  ab=tuple(ab);assert ab not in mp
  mp[ab]=sign*a
assert len(mp)==160
terms=[]
for (a,b),c in sorted(mp.items()):
 if a>b>=0 and c:
  terms.append((a,b,a*a+b*b,(2 if b==0 else 4)*c))
assert len(terms)==21
for (a,b),c in mp.items():
 assert mp[(-b,a)]==-c
 assert mp[(-a,-b)]==c
 assert mp[(a,-b)]==c
assert max(a+b for a,b,n,c in terms)==9
assert max(a for a,b,n,c in terms)==7
assert sum(abs(c) for a,b,n,c in terms)<37
assert set(n for _,_,n,_ in terms) == {1,4,5,9,10,13,16,17,20,25,26,29,34,36,37,40,41,45,49,50}
fromfile=[]
for line in (R/'terms.txt').read_text().splitlines():
 a,b,n,num,den=map(int,line.split());fromfile.append((a,b,n,Q(num,den)))
assert set(fromfile)==set(terms)
assert all(a*a+b*b==n for a,b,n,_ in terms)
# Tightened REAL Y radial interpolation (Lobatto degree 120) analytic envelope.
# rho=2, affine y=17/15+(8/15)z, real minimum on ellipse 7/15,
# |y|<=9/5, hence Re(2pi sqrt(n)y)>14/5 for n>=1.
# |K_ir(z)|<=exp(-a)*sqrt(pi/(2a))<1/16, likewise K'_ir,
# since a>14/5, exp(14/5)>16 and pi<22/7<2(a-1).
assert sum(Q(14,5)**k / __import__('math').factorial(k) for k in range(11))>16
assert Q(22,7) < 2*(Q(14,5)-1)
assert Q(22,7) < 2*Q(14,5)
assert Q(7)*Q(8)*Q(9,5) < 101  # 2pi sqrt(n) |y| <101
assert Q(9,5)/16 < Q(1,8)   # |F|<1/8
assert Q(1+101,16) < 7      # |Fy|<7
# Chebyshev-Lobatto uniform value/derivative error <=4 M rho^-120/(rho-1)
assert 4*Q(7)/2**120 < Q(1,10**32)
# Independent rational envelope for J on the 2D rho=3/2 polyellipse.
# Complex Chebyshev coordinates: x=s/2, t=(1+u)/4.
Uellipse=Q(13,24)**2+Q(25,48)**2
assert 1-Uellipse>Q(13,20)**2
assert 1+Uellipse<Q(13,10)**2
# Each Fourier product has |Im phase|<11 from 2pi<44/7,
# |Im x|<=5/24, |Im t|<=5/48, a<=7, a+b<=9.
assert Q(44,7)*Q(5,48)*(7+9)<11
from math import factorial
exp11=sum(Q(11)**k/factorial(k) for k in range(101))+Q(11)**101/factorial(101)/(1-Q(11,102))
assert exp11<60000
# Re(z)>39/10 and exp(39/10)>40: |K_ir|,|K'_ir|<1/40.
assert sum(Q(39,10)**k/factorial(k) for k in range(7))>40
# For F=yK and Fy=K+zK' use |y|<13/10, 2pi sqrt(n)<51.
assert Q(13,10)/40==Q(13,400)
assert Q(1,40)+Q(51)*Q(13,10)/40<Q(17,10)
# |Tx|,|Tt|<=2pi(a+b) exp(11)<7*9*60000.
Mcomplex=Q(13,10)*Q(37)*(Q(11,10)*Q(13,400)*Q(3780000)+Q(13,10)*Q(17,10)*Q(120000))
assert Mcomplex<Q(20_000_000)
# propagated J boundary error by < 37*(2pi*(a+b)+2)*1e-32,
# using 2pi<7, a+b<=9 and |T|<=2.
DELTA_J=Q(37)*(Q(63)+Q(2))*Q(1,10**32)
assert DELTA_J < Q(3,10**29)
# Interpolation on 128^2 boundary nodes, Gaussian Chebyshev tensor DCT error
N=128
assert 4*N*N*Q(3,10**29) < Q(2,10**24)
# Existing Bernstein ellipse analytic bound for actual J (the earlier
# independent mathematical audit remains an explicit hypothesis).
TAIL=144*Q(20_000_000)*Q(2,3)**128
assert TAIL<Q(8313,10**17) # 8.313e-14
# DCT coefficient intervals from MPFR outward verifier, use an exact
# rational ceiling with comfortable room above printed directed interval.
COEFF=Q(305402,10**17) # 3.05402e-12
assert COEFF==Q('3.05402e-12')
JUB=COEFF+Q(2,10**24)+TAIL
assert JUB<Q(314,10**14) # 3.14e-12
# Boundary hyperbolic area <5/6; trace constant 81/4.
LAMBDA=1+Q(6622119348586011,10**15)**2
assert LAMBDA<Q(45)
assert 1+LAMBDA+Q(8,10**6)<49
# Sign-only spectral test: if no eigenvalue is within eps then
# eps||f||<=sqrt(1+lambda+eps)* ||weak_residual|| <=7*(9/2)*JUB
# using ||f||>13e-6 and area_H<1.
EPS=Q(8,10**6)
assert 7*Q(9,2)*JUB < Q(13,10**6)*EPS
# MPFR table extrema, parsed as exact decimals; no claims about their
# provenance are derived merely from this parse.
getcontext().prec=100
widths={}
with (R/'radial_coeffs_192.csv').open() as fh:
 for row in csv.reader(fh):
  if row[0].startswith('#'):continue
  n=int(row[0]);a=Decimal(row[3])-Decimal(row[2]);b=Decimal(row[5])-Decimal(row[4]);
  w=widths.setdefault(n,[Decimal(0),Decimal(0)]);w[0]+=a;w[1]+=b
assert len(widths)==24
widthF=max(x[0] for x in widths.values());widthG=max(x[1] for x in widths.values())
assert widthF<Decimal('4e-36') and widthG<Decimal('2e-34')
assert Q(3000)*(Q(10,13)**384)<Q(1,10**38)
print('PASS: exact frozen Fourier projection, 21-term field and C input coefficients')
print('PASS: stronger radial rho=2 comparison and 192-point GL error inequality')
print('PASS: DCT interpolation and weak-residual arithmetic implication')
print('max radial coefficient interval widths F,Fprime:',widthF,widthG)
print('J analytic radial approximation error bound (before DCT):',float(DELTA_J))
print('J uniform ceiling from MPFR DCT plus interpolation:',float(JUB))
print('lambda exact decimal:',str(float(LAMBDA)))
print('epsilon spectral inclusion ceiling:',float(EPS))
print('frozen projected JSON SHA256:',hashlib.sha256(raw).hexdigest())
print('IMPORTANT: exact checks do not independently derive Picard geometry or entire analytic complex-ellipse lemmas')
