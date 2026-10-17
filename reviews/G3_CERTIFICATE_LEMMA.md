## G3. Analytic remainder lemmas (ρ_x = ρ_t = 3/2)

**Notation.** E_ρ is the Bernstein ellipse of parameter ρ about the relevant real interval. The tensor polyellipse is E_{3/2} × E_{3/2}, which is contained in |x| ≤ 13/24, |t| ≤ 25/48, |Im x| ≤ 5/24, |Im t| ≤ 5/48. Let w = 1 − x² − t² and y0 = √w.

**Lemma G3.1 (branch of y0).** On the closed polyellipse, |x² + t²| ≤ 1301/2304 < 0.565, so Re w ≥ 1 − 1301/2304 > 0.435. Hence w lies in the open right half-plane, and the principal square root y0 = √w is analytic on a neighbourhood of the closed polyellipse. It never meets w = 0 or the cut (−∞, 0]. Moreover, for Re w > 0 one has Re √w = √((|w|+Re w)/2) ≥ √(Re w) > √0.435 > 13/20, and |y0| = |w|^{1/2} < 13/10.

**Lemma G3.2 (complex Bessel bounds).** For real r and Re z > 0, K_ir(z) = ∫_0^∞ e^{−z cosh u} cos(ru) du and K'_ir(z) = −∫_0^∞ e^{−z cosh u} cosh u cos(ru) du. Hence |K_ir(z)| ≤ ∫_0^∞ e^{−Re z cosh u} du = K_0(Re z) and |K'_ir(z)| ≤ ∫_0^∞ e^{−Re z cosh u} cosh u du = K_1(Re z). Both bounds are decreasing in Re z. For Re z > 3.9 (here Re z ≥ 2π·(13/20)), Arb gives K_0(3.9) = 0.01248... and K_1(3.9) = 0.01400..., both < 1/40. (For the radial ρ = 2 region, K_0(2.8) = 0.0438 and K_1(2.8) = 0.0511, both < 1/16.)

**Lemma G3.3 (bound on J).** Using |F| ≤ 13/400, |F'| < 17/10, the largest phase factor over the 21 terms (e^{9.84} < e^{11} < 60000), and the triangle inequality, exact rational arithmetic gives sup|J| ≤ 1.870×10^7 < 2×10^7 on the polyellipse (script: `g3/exact_checks_output.txt`).

**Lemma G3.4 (tensor aliasing, N = 128 nodes per variable, degree 127).** If f is analytic and |f| ≤ M on E_ρ × E_ρ, then its Chebyshev coefficients satisfy |a_kl| ≤ 4M ρ^{−(k+l)}. Each index outside the 128×128 box contributes once as truncation tail and once as alias (first-kind nodes). With S = Σ_{k≥0} ρ^{−k} = 3 at ρ = 3/2, the error is at most 2·(4M)·2·S²·ρ^{−128} = 144·M·ρ^{−128}.

**Corollary.** With M = 2×10^7 and ρ = 3/2: 144 · 2×10^7 · (2/3)^128 = 8.3121×10^{−14} < 8.313×10^{−14} (exact rational check).

**Status.**
- Rigorous: G3.1, G3.2, the arithmetic in G3.3 and the Corollary.
- Cited: the coefficient bound |a_kl| ≤ 4Mρ^{−(k+l)}. The first-kind aliasing bookkeeping in G3.4 is our own derivation.
- Optional model-function test, still open: check numerically that the 144·M·ρ^{−N} bound holds for an explicit analytic function on the tensor grid.
- Optional: a validated ball-arithmetic sup of |J| on the torus. Arb's `bessel_k` over-enclosures imaginary-order ball arguments, so this needs a custom enclosure. It is not needed for the certificate (the sampled sup is about 3.5×10^{−12}).
