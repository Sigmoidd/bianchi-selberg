# Hybrid / periodize δ_aut vs M

r=6.0, Y0=0.8, jump_weight=2.0, n_face=16

| M | mode | δ_aut | τ_proxy | η_proxy | κ_eq | fac vs col | s |
|--:|------|------:|--------:|--------:|-----:|-----------:|--:|
| 48 | collocation | 3.386e+00 | 2.388e-01 | 3.625e+00 | 25.5 | 1.00 | 0.3 |
| 48 | multi | 4.431e-01 | 1.128e-01 | 5.559e-01 | 53.9 | 7.64 | 0.7 |
| 48 | periodize | 8.839e-02 | 1.128e-01 | 2.012e-01 | 53.9 | 38.31 | 1.5 |
| 64 | collocation | 1.621e+01 | 2.104e-01 | 1.642e+01 | 28.9 | 1.00 | 0.4 |
| 64 | multi | 8.620e-01 | 1.240e-01 | 9.861e-01 | 49.0 | 18.80 | 0.9 |
| 64 | periodize | 1.737e-01 | 1.240e-01 | 2.977e-01 | 49.0 | 93.34 | 2.5 |
| 100 | collocation | 2.587e+01 | 1.445e-01 | 2.601e+01 | 42.1 | 1.00 | 1.4 |
| 100 | multi | 1.713e+00 | 8.074e-02 | 1.794e+00 | 75.3 | 15.10 | 2.6 |
| 100 | periodize | 3.438e-01 | 8.074e-02 | 4.246e-01 | 75.3 | 75.23 | 7.6 |

Non-certifying. Hard map unchanged.
fac vs col = δ_collocation / δ_mode at same M.

## Honest note on larger M

Absolute δ_aut **increases** with M for this hybrid pin (best periodize still
M=48: δ≈0.088). Factors vs collocation stay large (10–90×), but bigger Fourier
truncation alone does not produce a near-automorphic trial. Next: true σ₀ /
conforming periodization, not M↑ on the present operator.
