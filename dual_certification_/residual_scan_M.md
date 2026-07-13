# Residual scan vs M

r=6.0, Y0=0.8, theta=0.5

| M | n_pts | n_modes | rel | τ_svd | κ_eq | η_proxy | s |
|--:|------:|--------:|----:|------:|-----:|--------:|--:|
| 100 | 484 | 316 | 2.375e-02 | 6.347e-01 | 42.1 | 1.682e-01 | 0.9 |
| 200 | 900 | 632 | 1.591e-02 | 5.970e-01 | 62.9 | 1.127e-01 | 2.7 |
| 400 | 1764 | 1256 | 1.032e-02 | 5.642e-01 | 96.9 | 7.309e-02 | 18.5 |
| 800 | 3364 | 2520 | 6.459e-03 | 5.092e-01 | 154.8 | 4.575e-02 | 137.9 |

log-log slopes: b(rel)≈-0.626, b(τ_disc)≈-42.984

Model collocation residual (single-cusp Hejhal-like). τ_disc = ||V a||₂ with unit-ℓ² near-kernel; rel = σ_min/σ_max after amp+Sinkhorn equilibration. Not Maass residual.

Non-certifying. Hard map unchanged.
