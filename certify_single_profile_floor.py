#!/usr/bin/env python3
"""Explicit interval attempt: Picard level-one Q-odd ONE Fourier profile.

Certified shape: analytic Gauss-Legendre error plus intervals for nodes,
weights, and elementary integral representation of K_{ir}, K'_{ir}.

MPMath IV is a prototype interval backend, NOT independently audited Arb.
All endpoint root tests are checked and an explicit error margin is added.
"""
import argparse, json, math, time
from pathlib import Path
import mpmath as mp
from scipy.special import roots_legendre

IV = mp.iv
R = mp.mpf('6.633')


def legendre_mpoly(n, x):
    p0, p1 = mp.mpf(1), x
    if n == 0: return p0, mp.mpf(0), p0
    for k in range(1, n):
        p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
    dp = n*(x*p1-p0)/(x*x-1)
    return p1, dp, p0

def legendre_ipoly(n, x):
    p0, p1 = IV.mpf(1), x
    if n == 0: return p0, IV.mpf(0)
    for k in range(1, n):
        p0, p1 = p1, ((2*k+1)*x*p1-k*p0)/(k+1)
    dp = n*(x*p1-p0)/(x*x-1)
    return p1, dp

def endpoint_str(x):
    return mp.nstr(x, 67, strip_zeros=False)

def nodes_and_weights(n, center='0', radius='1'):
    # Compute ordinary high-precision roots then *validate* each rational interval.
    centers, _weights = roots_legendre(n)
    rr = []
    roots = []
    eps = mp.mpf('1e-45')
    c=IV.mpf(center);rad=IV.mpf(radius)
    for d in centers:
        x=mp.mpf(str(d))
        for i in range(8):
            p,dp,_ = legendre_mpoly(n,x)
            x-=p/dp
        lo=endpoint_str(x-eps);hi=endpoint_str(x+eps)
        left = IV.mpf(lo);right=IV.mpf(hi)
        pl,_ = legendre_ipoly(n,left)
        pr,_ = legendre_ipoly(n,right)
        if not ((pl>0 and pr<0) or (pl<0 and pr>0)):
            raise ArithmeticError(f'Legendre root n={n} not bracketed: x={x} pl={pl} pr={pr}')
        box=IV.mpf([lo,hi]);_p,deriv=legendre_ipoly(n,box)
        weight=IV.mpf(2)/((IV.mpf(1)-box*box)*deriv*deriv)
        if not weight>0: raise ArithmeticError('negative/quasi-infinite Gauss weight')
        roots.append(box)
        rr.append((c+rad*box, rad*weight))
    for left,right in zip(roots[:-1], roots[1:]):
        if not left.b < right.a:
            raise ArithmeticError('Legendre root intervals are not disjoint and ordered')
    # Exact GL weight identity sum w_j=2, independently enclosed.
    total=sum((weight for (_x,weight) in rr), IV.mpf(0))
    if not (total.a <= 2*rad <= total.b):
        raise ArithmeticError('Gauss weights fail total-weight identity')
    return rr

def upper(x):return float(x.b)
def lower(x):return float(x.a)

def bessel_gl_series(u_nodes, y, r, epsilon):
    # K_{ir}(z) = integral_0^inf e^{-z cosh u} cos(r u) du
    # K'_{ir}(z) = -integral_0^inf e^{-z cosh u}cosh(u)cos(r u)du.
    z=IV.mpf(2)*IV.pi*y
    # Precomputed trig on validated u Gauss intervals
    val=IV.mpf(0); der=IV.mpf(0)
    for cu,coef,wt in u_nodes:
        q=wt*coef*IV.exp(-z*cu)
        val+=q
        der-=q*cu
    err=IV.mpf(['-'+str(epsilon),str(epsilon)])
    return val+err,der+err

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--space-n',type=int,default=56)
    ap.add_argument('--bessel-n',type=int,default=80)
    ap.add_argument('--bits',type=int,default=320)
    ap.add_argument('--out',default=str(Path(__file__).resolve().parent / 'single_profile_floor_interval_result.json'))
    args=ap.parse_args()
    mp.mp.dps=100
    IV.dps=max(95, math.ceil(args.bits/3.32192809488736))
    n=args.space_n;nu=args.bessel_n
    start=time.monotonic()
    # The Bessel integral Gauss error, using ellipse rho_u=13/10 and T=6,
    # is bounded by (80*M_u)*rho_u**(-(2*nu-1)), M_u<=2000.
    # Tail from 6 to infinity is <1e-100 (z>=4).
    u_error=IV.mpf('160000')*(IV.mpf(10)/13)**(2*nu-1)+IV.mpf('1e-100')
    # Use outward rounded decimal upper bound, derived independent of result.
    # Conservative constant 1e-12 exceeds the n_u=80 envelope 1.23e-13.
    if nu<80: raise ValueError('prototype bound requires bessel_n>=80')
    berr='1e-12'
    if not (u_error<IV.mpf(berr)):
        raise ArithmeticError(f'Bessel Gauss error exceeds selected eps: {u_error}')
    u_g=nodes_and_weights(nu,center='3',radius='3')
    r=IV.mpf('6.633')
    u_nodes=[(((IV.exp(u)+IV.exp(-u))/2),IV.cos(r*u),w) for u,w in u_g]
    print('validated Bessel Gauss nodes',nu,'time',round(time.monotonic()-start,1),flush=True)
    sx=nodes_and_weights(n,center='0',radius='0.5')
    st=nodes_and_weights(n,center='0.25',radius='0.25')
    print('validated spatial Gauss nodes',n,'time',round(time.monotonic()-start,1),flush=True)
    s=IV.mpf(0)
    for ix,(x,wx) in enumerate(sx):
        for t,wt in st:
            y=IV.sqrt(IV.mpf(1)-x*x-t*t)
            kval,kder=bessel_gl_series(u_nodes,y,r,berr)
            X=2*IV.pi*x;T=2*IV.pi*t
            U=IV.cos(X)-IV.cos(T)
            DU=2*IV.pi*(-x*IV.sin(X)+t*IV.sin(T))
            z=2*IV.pi*y
            # f=2yK U, - f (x f_x+t f_t+y f_y)/y²:
            integrand=-4*(kval*kval*(U*U+U*DU)+z*kval*kder*U*U)
            s+=wx*wt*integrand
        if ix%8==0: print('spatial row',ix,'/',n,'running B upper',upper(s),'sec',round(time.monotonic()-start,1),flush=True)
    # On rho=7/5 Bernstein ellipses in both spatial variables, the
    # analytic integrand has abs <=100: Re(y)>1/2, |y|<2,
    # |K|,|K'|<1/20, |U|<4, |DU|<15.
    # Tensor Gauss remainder <= 4M/(rho-1)*rho**(-(2n-1))
    # We widen by 10x for the independent 2D/product estimates.
    rho=IV.mpf(7)/5
    floor_error=IV.mpf(10000)*rho**(-(2*n-1))
    if not floor_error < IV.mpf('1e-12'):
        raise ArithmeticError('floor analytic bound exceeds 1e-12')
    out_interval=s+IV.mpf(['-1e-12','1e-12'])
    result={
        'status':'INTERVAL ATTEMPT, mpmath.iv backend (independent Arb replay required)',
        'r':'6633/1000', 'lambda0':'44996689/1000000',
        'gauss_space_n':n,'gauss_bessel_n':nu,'precision_bits':args.bits,
        'validated_legendre_root_brackets':2*n+nu,
        'bessel_quadrature_and_tail_error_upper':str(u_error),
        'per_K_or_Kprime_error_widening':berr,
        'analytic_floor_quadrature_remainder':str(floor_error),
        'quadrature_sum_enclosure':str(s),
        'floor_flux_enclosure':str(out_interval),
        'floor_flux_upper':str(out_interval.b),
        'sign_negative':bool(out_interval<0),
        'seconds':round(time.monotonic()-start,1),
    }
    with open(args.out,'w') as f:json.dump(result,f,indent=2)
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__': main()
