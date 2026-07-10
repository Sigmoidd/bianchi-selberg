"""Sensitivity scan: B(k, delta) across admissible test functions.
Consistency requirements: B >= 0 always (bounds a non-negative spectral sum);
B smooth in delta; conclusion B < 1 stable."""
import picard_stf as P
from flint import arb

L0 = float(((3 + arb(5).sqrt())/2).log())   # 0.96242...

print(f"systole l0 = {L0:.6f};  need 2k*delta <= l0")
print(f"{'k':>3} {'delta':>8} {'2k*delta':>9} {'B (certified enclosure)':>40}")
results = []
for k in (2, 3, 4):
    for frac in (0.999, 0.9, 0.75, 0.6):
        delta = L0 * frac / (2*k)
        B = P.compute_B(k, delta, R=80, prec=192, verbose=False)
        mid = (float(B.lower()) + float(B.upper()))/2
        print(f"{k:>3} {delta:>8.5f} {2*k*delta:>9.5f}   [{float(B.lower()):+.6f}, {float(B.upper()):+.6f}]  mid={mid:+.5f}")
        results.append((k, delta, float(B.lower()), float(B.upper())))

best = min(results, key=lambda t: t[3])
print(f"\nbest upper bound: B <= {best[3]:.6f} at k={best[0]}, delta={best[1]:.5f}")
neg = [r for r in results if r[3] < 0]
print("consistency (no negative B):", "FAIL" if neg else "PASS")
print("Selberg level 1 certified (all B<1):", "PASS" if all(r[3] < 1 for r in results) else "PARTIAL")
