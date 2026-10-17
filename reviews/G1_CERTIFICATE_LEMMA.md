## G1. Geometry and form domain (insert into FINAL_PICARD_BOUNDARY_CERTIFICATE.md)

**Setting.** Γ = PSL₂(Z[i]) acts on H³ = {(z,y): z = x+it, y>0}. Let
F = {−1/2 ≤ x ≤ 1/2, 0 ≤ t ≤ 1/2, x²+t²+y² ≥ 1}. This is the Picard Ford domain for Γ, folded by the reflection σ:(x,t)→(x,−t) (EGM Ch. 7 §7.3; Swan 1971). Its volume is 0.3053218647… (checked against the covolume). Let Q:(x,t)→(−t,x) be the quarter-turn. Q is **not** in Γ; it lies in the normalizer.

**Lemma G1.1 (pairing S).** Let S be the action of z ↦ −1/z. Then S(z,y) = (−z̄, y)/(|z|²+y²). It is a rotation by π about the geodesic {x=0, t²+y²=1}, and it is orientation-preserving. It is not the sector reflection σ. Restricted to the unit sphere |z|²+y²=1 it reads S(x,t,y) = (−x,t,y). It pairs the floor with itself, and J(Sp) = J(p).

**Lemma G1.2 (form domain).** The trial f is in H¹ of the Q-odd, reflection-even sector. Its Γ-extension is well defined because the traces of f agree across the pairings: x=±1/2 (translation), t=0 and t=1/2 (reflection-even, with integer frequencies giving f_t=0 there), and the floor (S). The proof checks the Fourier-Bessel terms on each face. The expansion decays like e^{−2πy} at the cusp. Stabilizers and edge-cycle angle sums were checked: order 4 at (0,0,1); 6 over z=1/2 and z=i/2; 12 at the corner (1/2,1/2,1/√2); 2 on the vertical axes and the fold x=0; 3 on the floor edges at x=±1/2 and t=1/2.

**Lemma G1.3 (weak residual identity).** For every v in the form domain,
a(f,v) − λ(f,v) = ∫_D J v̄ dx dt / y₀³, with y₀ = √(1−x²−t²), D = {|x|≤1/2, 0≤t≤1/2},
where J is the package's floor jump term. Proof outline: apply Green's identity on F, using (L−λ₀)f = 0 in the interior (exact).
- Vertical faces x=±1/2 and t=0,1/2: boundary terms vanish exactly for every v, since f_x=0 at x=±1/2 and f_t=0 at t=0,1/2.
- Cusp: the boundary term at height Y_k tends to 0 along a sequence Y_k→∞, with no assumption on the constant term of v.
- Floor: the outward normal is the radial direction toward the origin, with measure dx dt/y₀³. S pairs the floor with itself and its fluxes add, because J(Sp)=J(p). So the identity uses the full D, with the same full-F normalization for ‖f‖ and a.
Numerical corroboration (not part of the proof): (L−λ₀)f ≈ 3e-34 at a random interior point; a Green test with a non-symmetric v agrees to 12 digits; max J on a grid is 2.11e-12.

**Assumptions.** (A1) The Picard fundamental-domain theorem is taken from the literature. (A2) The form norm is ‖v‖² = a(v,v) + ‖v‖². (A3) The stabilizer enumeration searched only Gaussian-integer entries up to size 2, and the group names rely on EGM's classification.

**Citation items needing human verification** (not yet checked; these were written from memory):
- EGM Ch. 7 §7.3: exact theorem number and page numbers for the Picard domain. Only the section titles (Ch. 7, §8.8, §9.6) were verified.
- Swan 1971 (Advances in Math. 6): page numbers.
- Picard 1884, Bianchi 1892, Humbert 1919: theorem and page numbers.
- Grunewald–Schwermer: theorem and page numbers.
