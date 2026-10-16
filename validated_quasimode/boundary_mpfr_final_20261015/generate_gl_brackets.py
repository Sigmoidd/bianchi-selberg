"""Generate rational decimal brackets only; the C/MPFR verifier proves containment."""
from pathlib import Path
import sys
import mpmath as mp
from scipy.special import roots_legendre
mp.mp.dps=210
out=Path(__file__).parent

def leg(n,x):
 p0,p1=mp.mpf(1),x
 for k in range(1,n):p0,p1=p1,((2*k+1)*x*p1-k*p0)/(k+1)
 return p1,n*(x*p1-p0)/(x*x-1)

for n in map(int,sys.argv[1:] or [56,80]):
 x0,_=roots_legendre(n)
 roots=[]
 for d in x0:
  x=mp.mpf(str(d))
  for _ in range(12):
   p,dp=leg(n,x);x-=p/dp
  e=mp.mpf('1e-110') if n>=192 else mp.mpf('1e-40') if n>=88 else mp.mpf('1e-47')
  lo=mp.nstr(x-e,160);hi=mp.nstr(x+e,160)
  roots.append((lo,hi))
 path=out/f'gl_n{n}.txt'
 path.write_text(str(n)+'\n'+''.join(f'{l} {h}\n' for l,h in roots))
 print(path)
