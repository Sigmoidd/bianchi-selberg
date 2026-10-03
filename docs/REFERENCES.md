# References and their roles

Private OCR extracts are not runtime inputs. The following links and the
local notes identify the reproducible mathematical inputs without bundling
third-party source texts.

| Source | Stable reference | Role |
|---|---|---|
| Joshua S. Friedman, *The Selberg Trace Formula and Selberg Zeta-Function for Cofinite Kleinian Groups with Finite Dimensional Unitary Representations* (2005 thesis; 2006 preprint) | [arXiv:math/0612807v1](https://arxiv.org/abs/math/0612807v1) | Theorem 4.1.1, printed pp.41–42; full centralizer and E(R) definition on p.42. Cuspidal formulas and identity in §§4.3–4.4; spectral decomposition Theorem 3.8.1. |
| Balkanova, Chatzakos, Cherubini, Frolenkov, Laaksonen, *Prime Geodesic Theorem in the Three-dimensional Hyperbolic Space* | [arXiv:1712.00880v2](https://arxiv.org/abs/1712.00880v2) | Second trace-formula reference, Theorem 2.2; cross-check of conventions. |
| Elstrodt, Grunewald, Mennicke, *Groups Acting on Hyperbolic Space: Harmonic Analysis and Number Theory* (1998) | [Springer DOI:10.1007/978-3-662-03626-6](https://doi.org/10.1007/978-3-662-03626-6) | Historical classification citation, Ch.4 §4.3; centralizer/orbital context. The cited class count is retained as a legacy input, not newly independently verified from the book in M0. |
| Elstrodt, Grunewald, Mennicke, Astérisque 94 (1982), pp.43 onward | [Numdam original article](https://www.numdam.org/item/AST_1982__94__43_0/) | Related historical source. It is not a substitute for the 1998 classification citation. |
| Aurich, Steiner, Then, *Numerical computation of Maass waveforms and an application to cosmology* | [arXiv:gr-qc/0404020](https://arxiv.org/abs/gr-qc/0404020) | Picard numerical spectrum and quoted Matthies coefficients (eqs.83–86), used in `verify_matthies.py`; a numerical spectrum is not an interval proof. |
| Mehmet Haluk Şengün, *Arithmetic Aspects of Bianchi Groups* | [arXiv:1204.6697](https://arxiv.org/abs/1204.6697) | Arithmetic overview, volume and cusp context. |
| Alexander D. Rahm, *The homological torsion of PSL2 of the imaginary quadratic integers* | [arXiv:1108.4608](https://arxiv.org/abs/1108.4608) | Potential independent torsion cross-check; not used as new-field completeness evidence. |
| NIST DLMF, §5.11 | [Digamma asymptotics and error bounds](https://dlmf.nist.gov/5.11) | Context for the historical tail derivation; the current elementary series proof in `ANALYTIC_DERIVATIONS.md` avoids relying on a complex remainder factor. |

`SYSTOLES.md` supplies the trace minimization argument.
`ANALYTIC_DERIVATIONS.md` supplies the B-spline derivative, field constants,
contour reduction, prime-power normalization and tail bounds.
`NORMALIZATION_ISSUE.md` identifies exactly what the flip witness proves.
`D2_INVENTORY_PROOF.md` supplies the self-contained d=2 arithmetic proof and
the involution orbital derivation. `INVENTORY_PROOF.md` records the same
obligations still required for d=11,19. `D7_INVENTORY_PROOF.md` supplies
the complete d=7 argument. Neither new-field proof uses an external class-count table.

`D2_REVIEW_GUIDE.md` separates the general mathematical inputs, local
completeness arguments, finite replay checks and analytic enclosure.
`ADDING_GROUPS.md` gives the uniform implementation contracts.
