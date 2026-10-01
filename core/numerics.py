"""Shared precision and certified panel quadrature."""
from flint import arb, acb, ctx

PREC = 100
RTOL = 1e-5

def setp():
    ctx.prec = PREC

def integ(f,a,b,panel=4.0):
    a=arb(a); b=arb(b); n=max(1,int(float((b-a)/panel)+0.5)); step=(b-a)/n
    tot=acb(0)
    for i in range(n):
        tot += acb.integral(lambda z,_: f(z), a+i*step, a+(i+1)*step, rel_tol=RTOL)
    return tot
