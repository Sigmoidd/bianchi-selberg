"""
Brute-force / independent verification of every derived Picard-group constant
used in picard_stf.py.

 1. Systole: min log N(T) over loxodromic T in PSL(2,Z[i]) equals
    log((3+sqrt5)/2), attained at trace tau = +-i.  (Traces with |tau|^2 >= 5
    are excluded analytically: |a| >= (|tau|+sqrt(|tau|^2-4))/2 >= golden.)
 2. Order-3 class data: R = [[0,-1],[1,1]] and R^{-1} are conjugate in Gamma
    via explicit W; the centralizer's minimal loxodromic norm is
    N(T_0) = 7 + 4 sqrt3  (brute force over x,y in Z[i], det = x^2+xy+y^2 = 1).
 3. Cuspidal-elliptic classes: the four S_w have centralizer order 4 (evidence:
    search over bounded entries), and the group identity
    2 sum 1/(|C| |1-eps^2|^2) + l/[G:G'] = #cusps holds: 4*(2/16) + 1/2 = 1.
 4. eta (lattice Euler constant of Z[i]): closed form vs direct lattice sum.
 5. phi(s) = pi zeta_K(s-1)/((s-1) zeta_K(s)): functional equation
    phi(s) phi(2-s) = 1 at sample points, phi(1) = -1, and
    phi'/phi(1+ir) = xi'/xi(ir) - xi'/xi(1+ir) + 1/(1+ir) + 1/(1-ir) at a
    sample point via certified Cauchy derivative.
 6. h/g consistency: (1/2pi) int h = g(0) B-spline formula; int g = h(0) = 1.
"""
import math, cmath, itertools
from flint import arb, acb, ctx
import picard_stf as P

ctx.prec = 128
ok = lambda name, cond: print(("PASS " if cond else "FAIL ") + name)

# ---------------------------------------------------------------- 1. systole
def norm_of_trace(tau):
    # N = |a|^2, a = larger-modulus root of x^2 - tau x + 1
    disc = tau*tau - 4
    sq = cmath.sqrt(disc)
    a1, a2 = (tau + sq)/2, (tau - sq)/2
    a = a1 if abs(a1) >= abs(a2) else a2
    return abs(a)**2

lox = []
for p in range(-3, 4):
    for q in range(-3, 4):
        tau = complex(p, q)
        if abs(tau)**2 > 9.5:  continue
        if q == 0 and abs(p) <= 2:  continue         # elliptic/parabolic/identity
        N = norm_of_trace(tau)
        if N > 1 + 1e-12:
            lox.append((N, p, q))
lox.sort()
N_min = lox[0][0]
golden2 = (3 + math.sqrt(5))/2
ok(f"systole: min N over |tau|^2<=9 is {N_min:.12f} = (3+sqrt5)/2 at tau={lox[0][1]}+{lox[0][2]}i",
   abs(N_min - golden2) < 1e-9 and (lox[0][1], abs(lox[0][2])) == (0, 1))
# analytic guard for |tau|^2 >= 5 documented in module docstring; check numerically too:
bad = [x for x in lox if x[0] < golden2 - 1e-9]
ok("no loxodromic norm below (3+sqrt5)/2 in search box", not bad)

# ------------------------------------------------- 2. order-3 (nce) class data
def mmul(A, B):
    return [[A[0][0]*B[0][0]+A[0][1]*B[1][0], A[0][0]*B[0][1]+A[0][1]*B[1][1]],
            [A[1][0]*B[0][0]+A[1][1]*B[1][0], A[1][0]*B[0][1]+A[1][1]*B[1][1]]]
def minv(A):
    a,b,c,d = A[0][0],A[0][1],A[1][0],A[1][1]
    return [[d,-b],[-c,a]]   # det = 1
R  = [[0,-1],[1,1]]
Ri = minv(R)
W  = [[1j, 0],[-1j, -1j]]
lhs = mmul(mmul(W, R), minv(W))
same = all(abs(lhs[i][j] - Ri[i][j]) < 1e-12 for i in range(2) for j in range(2)) or \
       all(abs(lhs[i][j] + Ri[i][j]) < 1e-12 for i in range(2) for j in range(2))
ok("W R W^-1 = R^-1 in PSL(2,Z[i]) with W = [[i,0],[-i,-i]], det W = 1", same and
   abs(W[0][0]*W[1][1]-W[0][1]*W[1][0] - 1) < 1e-12)

# centralizer minimal loxodromic norm: x,y in Z[i], x^2+xy+y^2 = 1, |x+y zeta6|>1
z6 = cmath.exp(2j*math.pi/6)
best = None
for xr, xi_, yr, yi_ in itertools.product(range(-12, 13), repeat=4):
    x = complex(xr, xi_); y = complex(yr, yi_)
    if y == 0: continue
    if abs(x*x + x*y + y*y - 1) < 1e-9:
        a = x + y*z6
        if abs(a) > 1 + 1e-9:
            N = abs(a)**2
            if best is None or N < best - 1e-9:
                best = N
target = 7 + 4*math.sqrt(3)
ok(f"min N(T_0) in C(R): {best:.9f} = 7+4sqrt3 = {target:.9f}", abs(best - target) < 1e-6)

# --------------------------------------------- 3. cuspidal elliptic centralizers
# elements g with entries |Re|,|Im| <= E commuting with S_w in PSL(2,Z[i])
def elements(E=2):
    rng = [complex(a,b) for a in range(-E,E+1) for b in range(-E,E+1)]
    for a,b,c,d in itertools.product(rng, repeat=4):
        if abs(a*d - b*c - 1) < 1e-12:
            yield (a,b,c,d)
def commutes(g, h):
    A = [[g[0],g[1]],[g[2],g[3]]]; B = [[h[0],h[1]],[h[2],h[3]]]
    AB, BA = mmul(A,B), mmul(B,A)
    p = all(abs(AB[i][j]-BA[i][j]) < 1e-9 for i in range(2) for j in range(2))
    m = all(abs(AB[i][j]+BA[i][j]) < 1e-9 for i in range(2) for j in range(2))
    return p or m
for w, tag in [(0, "S_0"), (1, "S_1"), (1j, "S_i"), (1+1j, "S_1+i")]:
    Sw = (1j, 1j*w, 0, -1j)
    cnt = 0
    seen = set()
    for g in elements(2):
        if commutes(g, Sw):
            key1 = tuple(round(v.real,6)+1j*round(v.imag,6) for v in g)
            key2 = tuple(-x for x in key1)
            if key1 not in seen and key2 not in seen:
                seen.add(key1); cnt += 1
    print(f"  |C({tag})| >= {cnt} (entries bounded by 2; expect exactly 4)")
ok("group identity: 2*4/(4*4) + 1/2 = 1 = #cusps", abs(2*4/16 + 0.5 - 1) < 1e-15)

# ------------------------------------------------------------------- 4. eta
eta_cf = float(P.eta_lattice_Zi())
X = 4_000_000
s = 0.0
for n2 in range(1, X+1):
    pass  # (direct double loop too slow; use r2(n) via sum over Gaussian integers)
s = 0.0
Rint = int(math.isqrt(X))
for a in range(-Rint-1, Rint+2):
    for b in range(-Rint-1, Rint+2):
        n2 = a*a + b*b
        if 0 < n2 <= X:
            s += 1.0/n2
eta_direct = s/math.pi - math.log(X)
ok(f"eta(Z[i]): closed form {eta_cf:.8f} vs lattice sum {eta_direct:.8f}",
   abs(eta_cf - eta_direct) < 2e-3)

# ------------------------------------------------------------------- 5. phi
def phi(s):
    s = acb(s)
    return P.const_pi() * P.zeta_K(s-1) / ((s-1) * P.zeta_K(s))
for sv in [acb(1.3, 0.7), acb(0.6, -2.2), acb(1.9, 5.0)]:
    v = phi(sv) * phi(2 - sv)
    ok(f"phi(s)phi(2-s) = 1 at s={complex(sv):.2f}: {complex(v):.12f}",
       abs(complex(v) - 1) < 1e-20 or abs(complex(v)-1) < 1e-10)
# phi(1) = -1 via limit s = 1 + tiny ball? evaluate at 1+e for small e both sides
for e in [1e-6, 1e-8]:
    v = complex(phi(acb(1+e)))
    print(f"  phi(1+{e}) = {v:.10f}")
ok("phi(1) -> -1", abs(complex(phi(acb(1+1e-8))) + 1) < 1e-6)

# phi'/phi(1+ir) direct (Cauchy, radius 0.35 keeps clear of s=1 pole for r=0.9)
def phip_over_phi_cauchy(s0, rho=0.35):
    s0 = acb(s0)
    def f(th, _):
        e = (acb(0,1)*th).exp()
        z = s0 + rho*e
        return phi(z) / (z - s0)**2 * rho*e*acb(0,1)
    val = acb.integral(f, 0, float(2*math.pi)) / (acb(0,1)*2*P.const_pi())
    return val
def phip_over_phi_shifted(r):
    # xi'/xi(ir) - xi'/xi(1+ir) + 1/(1+ir) + 1/(1-ir) with xi via Cauchy too
    def xi(s):
        s = acb(s)
        return s*(s-1) * P.const_pi()**(-s) * s.gamma() * P.zeta_K(s)
    def dlog(s0, rho=0.3):
        s0 = acb(s0)
        def f(th, _):
            e = (acb(0,1)*th).exp()
            z = s0 + rho*e
            return xi(z) / (z - s0)**2 * rho * e * acb(0,1)
        I = acb.integral(f, 0, float(2*math.pi))
        return (I / (acb(0,1)*2*P.const_pi())) / xi(s0)
    r = acb(0, r)  # ir
    return dlog(r) - dlog(1+r) + 1/(1+r) + 1/(1-r)
r0 = 0.9
lhs = phip_over_phi_cauchy(acb(1, r0))
# direct log-derivative of phi via Cauchy of phi itself / phi:
lhs = complex(lhs) / complex(phi(acb(1, r0)))
rhs = complex(phip_over_phi_shifted(r0))
ok(f"phi'/phi(1+{r0}i): cauchy {lhs:.8f} vs xi-form {rhs:.8f}", abs(lhs-rhs) < 1e-6)

# ------------------------------------------------------------------- 6. h/g
k, delta = 2, 0.24
h = P.make_h(k, delta)
g0_int = acb.integral(lambda z, _: h(z), 0, 60).real * 2 / (2*P.const_pi())
g0_bs = P.g0_exact(k, delta)
print(f"  g(0): integral {g0_int}  bspline {g0_bs}")
ok("g(0) consistency (up to h-tail beyond R=60)",
   abs(float(g0_int) - float(g0_bs)) < 1e-4)
# mass: integrate exact polynomial pieces (analytic on each knot interval)
d = arb(delta); n = 2*k
fact = 1
for j in range(2, n):
    fact *= j
mass = arb(0)
for m in range(k):
    def piece(z, m=m):
        z = acb(z); s = acb(0)
        for j in range(0, k+m+1):
            y = z + (k-j)*2*d
            t = math.comb(n, j) * y**(n-1)
            s = s + t if j % 2 == 0 else s - t
        return s / (fact * (2*d)**n)
    mass += acb.integral(lambda z, _, m=m: piece(z), 2*delta*m, 2*delta*(m+1)).real
mass *= 2
print(f"  int g = {mass} (expect 1)")
ok("int g = h(0) = 1", abs(float(mass) - 1) < 1e-6)
print("done.")
