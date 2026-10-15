# BH? No: Picard level-one narrow inclusion, boundary-only final verification

**Repo:** `Sigmoidd/bianchi-selberg`, draft PR #14. **Today's novel branch should be reviewed independently before merging.**

Objective: certify **one** odd Picard cuspidal eigenvalue in `(44.85246466691721-1e-5,44.85246466691721+1e-5)` without any three-dimensional Poincaré-average residual. The center is exactly `lambda=1+(6622119348586011/10^15)^2`; compute interval endpoints rationally rather than rounding text.

## Accepted mathematical reduction (requires independent referee check)

The *exact* frozen 21-term reflection-even/quarter-turn-odd Fourier–Bessel seed is H1-automorphic by **value matching** across all paired Ford faces but does not have matching normal derivatives. Its weak residual is the Ford-floor normal derivative `J=-y(x f_x+t f_t+y f_y)` on `y=sqrt(1-x²-t²)`. With a rational trace constant `C=81/4`, area `<5/6`, and previously MPFR-certified `||f||>13e-6`, proving `||J||_L∞ <4e-12` suffices for a `10^-5` spectral inclusion.

## Exact numerical gate

1. Regenerate frozen `f` from `inputs/analytic_quasimode_coefficients_50_reflection_even.json` and `inputs/radial_coeffs_all.csv` and independently verify SHA/normalization.
2. Construct **128×128 Gauss-Chebyshev nodes** on `x=s/2`, `t=(1+u)/4`. Use MPFR outward arithmetic for pi, cos/sin, the square root, the 24 norm radial polynomials and their first derivative. Enclose exact nodal J values using the certified radial interpolation error `<=1e-20` in F and `<=4e-20` in F'.
3. Compute 2D DCT-II (degree127×127) with MPFR outward intervals. Verify **sum of absolute coefficient suprema < 3.5e-12**. The current 128×128 ordinary IEEE diagnostic is `3.0540183168045984e-12`; NOT a certificate.
4. Re-check the analytic polyellipse bound `|J|<20e6` for `rho=3/2` using explicit complex Bessel integral inequalities. The double Chebyshev alias remainder is `<144*(20e6)*(2/3)^128 <1e-13` (audit_theorem.py checks rational part).
5. Conclude uniform J `<3.6e-12<4e-12`, apply proved trace/spectral theorem, output certified enclosure endpoints and machine-replayable SHA logs. Use a second independent verifier or mutual check of interval FFT normalization.
6. Test mutations: J sign, omit lowest Fourier term, perturb 1% coefficient, r changed by 1e-4, lost coefficient interval radius, erroneous H1 face conditions.

## Hard review conditions

- Independently verify side-face pairings, *reflection-even sector versus actual PSL2(Z[i])* and orbifold trace conventions.
- Do not confuse `J` weak boundary residual with the L2 residual of the Poincaré average `W`. Those are two different a posteriori methods.
- `d<eps` proves SOME eigenvalue, not uniqueness, index, multiplicity, and no analogue at Gamma0(2+i) or trace-formula completeness.
- Do not mark existing theorem flags certified until MPFR DCT and geometry have both passed independent review.
- Free/local compute only.

## Existing scripts and verified numbers

- `audit_theorem.py`: exact fractions for trace, Cauchy inclusion, polyellipse, Taylor exponential constants.
- `boundary_weak_residual.py`: reconstruction and independent double boundary Gauss + Chebyshev diagnostics.
- `independent_direct_bessel_spots.py`: direct mpmath Bessel checks at three distinct floor locations.
- `test_boundary.py`: five local tests.
- `BOUNDARY_METHOD_EXPERT_REVIEW.md`: complete mathematical derivation and reference/limitation audit.

**Stop criterion:** If rigorous coefficient l1 >=3.5e-12, refine to 160 nodes or certify boundary L2 by validated 2D quadrature. Never substitute a sampled maximum for global bound.
