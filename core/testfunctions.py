"""Certified sinc and piecewise-polynomial cuspidal integral."""
from flint import arb, acb
import math

def sinc2k(z,k,d):
    # Arb's entire sinc handles zero and complex balls directly; no midpoint
    # branch or custom Taylor-remainder contract is needed.
    return (acb(z)*arb(d)).sinc()**(2*k)
def h_isig(k,d,sig=1):
    return sinc2k(acb(0, arb(sig)), k, d).real

def ce_integral(k,d,supp,ckern):
    """int_0^supp g(x) sinh x/(cosh x+ckern) dx; polynomial per knot interval."""
    d=arb(d); n=2*k; fact=arb(1)
    for j in range(2,n): fact=fact*j
    tot=acb(0)
    for m in range(k):
        a_,b_=2*d*m,2*d*(m+1)
        def piece(z,m=m):
            z=acb(z); s=acb(0)
            for j in range(0,k+m+1):
                y=z+(k-j)*2*d
                t=math.comb(n,j)*y**(n-1)
                s = s+t if j%2==0 else s-t
            gp=s/(fact*(2*d)**n)
            return gp*z.sinh()/(z.cosh()+ckern)
        tot += acb.integral(lambda z,_: piece(z), a_, b_)
    return tot.real
