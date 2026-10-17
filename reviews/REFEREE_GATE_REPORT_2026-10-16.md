# Picard Q-odd boundary spectral inclusion: hostile-referee checkpoint

Date: 2026-10-16. Independent working session. **Verdict: arithmetic evidence reproduced; final theorem not independently certified.**

## Precise proposed claim

For Gamma = PSL_2(Z[i]) in the Q-odd, reflection-even cusp sector, put
r = 6622119348586011/10^15 and lambda0=1+r^2. The claim is existence of at least one
cuspidal eigenvalue in (lambda0-8e-6, lambda0+8e-6). This makes no claim of uniqueness,
firstness, or a congruence subgroup.

## Verified in this session

- Extracted the immutable 22-file October 15 package and verified all 20 SHA-256 manifest entries against file bytes: PASS.
- Ran `python audit_exact_theorem.py`: PASS, including frozen 21-term projection and rational error-budget implications.
- Ran `python audit_seed_slab_simpson_bounds.py`: PASS, including finite-interval remainder inequalities.
- Compiled the independently replayable MPFR C DCT verifier with `cc -O2 -std=c11 -Wall -Wextra` against libmpfr.so.6 and libgmp. Official MPFR development headers are absent, so code selected its bundled MPFR-ABI fallback. Compiler warnings were limited to misleading indentation and unused functions; no diagnostic asserted a result mismatch.
- Fresh 8-node run: outward coefficient L1 upper 3.018981521794834110078e-12, PASS.
- Fresh 32-node run: outward coefficient L1 upper 3.054010740670342054988e-12, PASS.
- Verified saved full 128-node output SHA through manifest. It records outward coefficient L1 upper 3.054010740798494611700e-12 and `CERTIFICATE_DCT_L1_PASS`. Did not rerun full 128-node computation in this independent review.
- Independently calculated exact-rational last step: 7*(9/2)*(3.14e-12) < (13e-6)*(8e-6), with ratio 0.9510577 and margin 5.09e-12. This is a **conditional** implication only.

## Referee blockers requiring a separate mathematical derivation

**G1. Automorphic H1 form domain and geometry.** Confirm from a presentation/fundamental domain that x=±1/2, t=0,1/2 (via the reflection-even sector), and Ford sphere S give exactly the asserted trace identities for this function; check stabilizers/corners and double-counting of the floor. Re-derive the equality `a(f,v)-lambda(f,v) = integral_floor J v` for every form-domain test v, including orientation and cusp limits.

**G2. Functional-analytic spectral lemma.** Prove in the actual Q-odd, reflection-even Hilbert space that the cusp constant terms vanish at all relevant cusps, the form embedding is compact, and the form-dual norm controls spectral distance with the displayed `(1+lambda_j)^(-1)` weights. The spectral inequality is elementary once this context is established.

**G3. Analytical complex-neighborhood estimates.** Verify independently the rho=3/2 *two-variable* Bernstein-ellipse bounds, including the branch of sqrt(1-x^2-t^2), the complex-order Bessel integral bound, and the tensor Chebyshev aliasing factor 144. The exact-rational script checks inequalities *given* these lemmas, not proofs of the lemmas.

**G4. Certified radial input provenance.** Cross-check that the MPFR radial quadrature and Chebyshev reconstruction enclose the exact 21-term seed and its radial derivatives, including all truncation and interpolation errors, and that the frozen data hashes match code inputs.

**G5. Truly independent arithmetic implementation.** Reproduce the 128x128 boundary node/DCT enclosure with a separately coded Arb/Acb or interval-validated Chebyshev evaluation, and build the original code against actual `mpfr.h` headers if available. The saved full run is reproducible in principle but the present session only freshly reran N=8 and 32.

## Scope and action

No valid basis to mark this publication-grade theorem or unique/first eigenvalue yet. Do not alter existing Rung-4 proof flags. Freeze optimization and use one PR checklist for G1-G5; if any gate fails, reopen that gate only. The most consequential mathematical risk is G1/G3, not extra decimal digits in the DCT.
