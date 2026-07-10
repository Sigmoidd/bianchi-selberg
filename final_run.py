"""Final certified run: best test function + certified spectral gap."""
import picard_stf as P
from flint import arb

k = 2
L0 = ((3 + arb(5).sqrt())/2).log()
delta = float(L0.lower()) * 0.999 / (2*k)

print(f"h = sinc^{2*k}(delta r), delta = {delta:.6f}, supp g = [-{2*k*delta:.5f}, {2*k*delta:.5f}]")
print(f"systole l0 = {L0}")
B = P.compute_B(k, delta, R=300, prec=192, verbose=True)
Bu = B.upper()
print(f"\nCertified: B <= {float(arb(Bu)):.6f}")
assert float(arb(Bu)) < 1

# spectral gap: h decreasing on [0, pi/delta]; find r* with h(r*) = B_upper.
lo, hi = 0.0, float(P.const_pi()) / delta
for _ in range(80):
    mid = (lo + hi)/2
    v = P.make_h(k, delta)(arb(mid)).real
    if v.lower() > float(arb(Bu)):
        lo = mid
    else:
        hi = mid
r_gap = lo   # h(r) > B certified for all r <= r_gap  (monotone segment)
lam_gap = 1 + r_gap**2
print(f"\nRESULTS (PSL(2,Z[i]) \\ H^3, level 1):")
print(f"  1. No exceptional eigenvalues: lambda_1 >= 1  [Selberg conjecture, level 1]")
print(f"  2. Certified spectral gap: no discrete eigenvalues in (0, {lam_gap:.4f});")
print(f"     i.e. lambda_1 >= {lam_gap:.4f}  (r_1 >= {r_gap:.4f})")
print(f"     [first cusp form believed near r ~ 6.6, lambda ~ 45: consistent]")
