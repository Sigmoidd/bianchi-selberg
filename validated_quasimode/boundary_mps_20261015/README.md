# Picard boundary-tension spectral inclusion research

This package is a **mathematical and numerical checkpoint, not a certified two-sided theorem**. It replaces the three-dimensional analytic residual of the automorphic Poincare-average `W` with the weak boundary residual of the exact interior Fourier-Bessel seed `f`.

Run on a standard free scientific Python environment:

```bash
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_boundary
python audit_theorem.py
OPENBLAS_NUM_THREADS=1 python boundary_weak_residual.py
python independent_direct_bessel_spots.py
```

The exact conditions, supporting arguments and essential caution flags are in `BOUNDARY_METHOD_EXPERT_REVIEW.md`. The precise next coding assignment is `NEXT_AGENT_PROMPT.md`.

**The remaining theorem gate is an outward-rounded MPFR DCT on 128 × 128 Chebyshev nodes**. Its coefficient absolute-sum must be <3.5e-12. Ordinary double evaluations consistently return ~3.054e-12 but do not certify the bound. The exact polyellipse tail (<8.32e-14) is independently audited in `audit_theorem.py` conditional on the complex Bessel inequalities in the review. When the nodal interval DCT and all face identifications have been independently checked, the proof may assert a Q-odd eigenvalue within 1e-5 of lambda=1+r², r fixed as 6.622119348586011. No eigenvalue index, uniqueness or trace-formula completeness follows.
