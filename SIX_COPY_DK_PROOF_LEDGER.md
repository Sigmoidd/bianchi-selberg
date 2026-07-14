# Six-copy / two-cusp D(K) proof ledger

**Status:** theorem architecture established; `rung4_certified=false`.

The central change is conceptual.  The six corrected components are not six
independent scalar functions.  They are the induced permutation local system

\[
L^2\!\left(\Gamma_0(2+i)\backslash\mathbb H^3\right)
\simeq
L^2\!\left(\mathrm{PSL}_2(\mathbb Z[i])\backslash\mathbb H^3;
\rho_{\rm ind}\right).
\]

This equivalence is exact, unitary, and Laplacian-intertwining.  It turns the
two subgroup cusps into two singular channels of one rank-six parabolic
holonomy.

## Dependency graph

```text
exact F5 right-coset action
        |
        v
unitary induced local system -- exact norm/Laplacian equivalence
        |
        +--> two parabolic orbits = two cusp channels
        |
        v
finite Whittaker trial (termwise PDE residual = 0)
        |
        v
certified C2 chart partition + active transitions
        |
        +--> value overlap defect delta0 in C^6
        +--> gradient overlap defect delta1 in C^6
        +--> partition constants b0, b1
        |
        v
exactly automorphic quasimode U
        |
        +--> target route: oldspace projector + Picard C-symmetry projector
        |                         |
        |                         v
        |                   purely cuspidal sector
        |
        +--> general route: normalized Hecke cuspidalizer (still theoretical)
                                  |
                                  v
                         purely cuspidal sector
        |
        v
spectral theorem: |lambda_j-lambda| <= residual / projected norm
```

## Proved or exactly verified (GREEN)

| Item | Evidence |
|---|---|
| Glue convention | `pi_g(c)=c.g^{-1}` agrees with the corrected operator identity |
| Induced representation | `(rho(g)v)_c=v_(c.g)` is a unitary permutation representation |
| Scalar/vector equivalence | union of six level-one tiles is a subgroup fundamental domain; norms match with no factor `1/6` |
| Laplacian equivalence | coset representatives act isometrically |
| Cusp count | parabolic and translation images have orbits `{0}` and `{1,2,3,4,5}` |
| Singular degree | exactly `2`; regular-channel dimension exactly `4` |
| Full finite image | exact order `60`, transitive on all six cosets |
| `TiR` role | already in `<T1,R>`; retained as a stored-transition cross-check |
| Symbolic PDE residual | every finite Whittaker term solves `(Delta-(1+r^2))F=0` |
| Smooth-gluing inequality | product rule gives `R <= tau + ||B0 d0 + B1 d1||_2` |
| Spectral enclosure after cuspidal projection | direct self-adjoint spectral expansion |
| Published target class | layout-preserving extraction of Then's table puts `r=6.62211934...` in class `C` |

Exact finite algebra is reproduced by `six_copy_DK_algebra.py` and
`six_copy_DK_algebra_result.json`.  The script uses integers and rational
numbers only.

## Evidence for the target-specific shortcut (YELLOW, strong)

Then's class `C` satisfies

\[
f(iz,y)=-f(z,y),\qquad f(-\bar z,y)=f(z,y).
\]

The quarter-turn fixes the level-one cusp datum, so the level-one Eisenstein
and residual spectrum is quarter-turn even.  Therefore the exact projector

\[
P_{\rm target}=
\frac{I-Q}{2}\,
\frac{\mathbf1_6\mathbf1_6^*}{6},
\qquad Q(z,y)=(iz,y),
\]

has cuspidal range.  This avoids the general Hecke-normalization problem for
the current mode.

There is no hidden simplicity hypothesis here. Although `Q:z -> iz` has
order four on hyperbolic space, `Q^2:z -> -z` is the action of
`diag(i,-i)` in `PSL2(Z[i])`. Hence the pullback on level-one automorphic
functions has square one, is self-adjoint, and `(I-Q)/2` is exactly
idempotent. The possible `+/- i` eigenspaces from an abstract order-four
operator do not occur on this quotient.

The fixed `M=8` coefficient vector independently supports this choice:

| Diagnostic | Relative error |
|---|---:|
| `a_(i beta) = -a_beta` | `1.365e-3` |
| `a_(conj beta) = a_beta` | `1.276e-3` |
| zero-cusp mass off the level-one lattice | `5.204e-4` |
| matched two-cusp coefficient mismatch | `6.076e-4` |
| sampled oldspace mass fraction (Q-invariant set) | `0.9999998893` |
| sampled oldspace/quarter-turn-odd mass fraction | `0.9999998893` |
| certified Track-B projected mass | `mu_B > 0.20380681000849094` on a cusp plateau |

These are floating diagnostics from `six_copy_target_symmetry.py`; they are
not proof inputs.  The missing proof input is an Arb lower bound for
`||P_target U||_2` after exact automorphization.

The exact target projector used by the enclosure is rotation odd, hence it
contains Then's `C` and `H` sign sectors. This is deliberate: either sector
is cuspidal, and the spectral conclusion does not assign a class label. An
additional commuting reflection-even projector would isolate class `C`, but
is not needed for certification and could only reduce the projected mass.

### Certified Track-B mass witness

`track_b_projected_mass_arb.py` collapses the five translated zero-cusp
copies exactly onto the integral infinity-cusp lattice, averages under the
Picard relation `z -> -z`, and then applies the rotation-odd projector. On
the cusp slab

`[-1/2,1/2] x [0,1/2] x [1.20,1.45]`

it proves

\[
\|W_B\|_2>0.20380681000849094.
\]

at 192-bit precision with 512 interval height segments. Full-torus Parseval
divided by two and a direct half-torus Rump Gram computation independently
prove positivity. A direct six-component evaluation agrees with the collapsed
Fourier formula in Arb; the largest difference radius is `5.07e-51`.

The exact quintic cusp cutoff is supported above `y=1.01` and equals one for
`y>=1.20`. Every other weight is `(1-chi_B) phi_j` for a normalized
subordinate family `sum phi_j=1`. Therefore the glued field equals `W_B` on
the witness slab and `epsilon_B=0` exactly. The
projected-mass input of Theorem D(K) is now closed. The current absolute
Taylor majorant from `continuum_defect_arb.py` remains irrelevant to this
mass proof; residual certification is still outstanding.

## Why the inherited one-cusp D(K) proof is not a safe base

Two steps in the current `theorem_DK.tex` need replacement, independently of
the number of cusps.

1. It starts with a scalar pointwise value defect and then invokes a reverse
   trace extension into `H1`, followed by `H2` control of its Laplacian.
   A bounded extension requires the correct fractional trace regularity (and
   normal compatibility for an `H2` patch); an `L2` value jump is not enough.

2. Above the continuous threshold, projecting a small-residual function onto
   the cuspidal spectrum can destroy all of its norm.  A continuous spectral
   wave packet can have arbitrarily small residual and no cusp-form component.
   Booker--Strombergsson--Venkatesh use an actual cuspidalizing operator; the
   soft statement in the inherited Step 5 does not supply the needed lower
   bound.

The new smooth partition lemma fixes the first issue by requiring both value
and first-gradient overlap defects.  The target symmetry projector fixes the
second issue for `r=6.622119...`.

## Independent review and repairs

An adversarial second pass found and repaired two substantive issues:

- Forward right actions are anti-homomorphic.  From `TiR = Ti R` one gets
  `f_Ti = f_R^{-1} o f_TiR`, so the exact verifier now reports
  `Ti = (0,4,5,1,2,3)`.  The previous value was its inverse.  The generated
  parabolic subgroup, its two orbits, singular degree two, and projector are
  unchanged and are rechecked from scratch by the script.
- The gluing lemma now compares candidates only after explicit unitary
  flat-bundle transport, uses covariant gradients, requires weighted local
  sections to extend `C2` by zero, restricts the norm anchor to an active
  local candidate, and proves membership in the self-adjoint Laplacian domain
  from `L2` distributional control.  These are proof hypotheses, not optional
  implementation details.

The review also confirmed that value and first-gradient overlap defects are
sufficient: partition identities cancel the common value and gradient terms,
so no independent second-derivative jump bound is required.

A second pass identified the elliptic-isotropy condition hidden in the phrase
"local section." The theorem now constructs one exactly by averaging the raw
Whittaker candidate over each finite chart stabilizer. Isotropy averaging
commutes with the Laplacian, and its value/gradient cost is included in the
same certified overlap defects. The target projector is now stated to commute
with the Laplacian, and the Hecke branch includes its explicit
projected-residual inequality.

## Remaining certificate inputs (YELLOW)

| Input | Why missing | Smallest provable addition |
|---|---|---|
| finite-type `C2` chart partition | no explicit weights or derivative enclosures exist | build a thickened Humbert cover with exact active transition list |
| `b0=||sum |Delta psi_j|||_2` | depends on that partition | Arb integrate polynomial/rational derivative bounds |
| `b1=||2 sum |grad psi_j|||_2` | same | same certificate |
| vector `delta0` | current code reports a scalar-component max | evaluate all six components and take the Arb Euclidean norm |
| vector `delta1` | current continuum code only bounds values | interval automatic differentiation of `F_c-rho(g)F(gP)` |
| edge/vertex transitions | four generators do not enumerate all chart overlaps | certify each active word directly, or a fully contained telescoping chain |
| elliptic uniformizers | raw approximate automorphy is not stabilizer equivariance | enumerate each finite isotropy group, average exactly, and include the averaged value/gradient defects |
| two-cusp domain/tail proof | zero modes and transported decay are not enclosed at both scaled cusps | certify zero constant terms, exponential `L2` decay, cusp-periodic local finiteness, and integrability of `B0 d0+B1 d1` |
| useful cancellation-preserving enclosures | first-order absolute modal sums give `5.05e4` | Taylor models or validated spectral interpolation on each active overlap |

The theorem output for the target route is

\[
|\lambda_j-\lambda|
\le
\frac{\tau+b_0\delta_0+b_1\delta_1}{\mu_{\rm target}},
\qquad \tau=0.
\]

There is no artificial `eta0`; the only smallness condition is the necessary
one `mu_target>0`.

## General two-cusp branch (RED until normalized)

For modes not isolated by a symmetry sector, the draft introduces

\[
\Diamond_{\mathfrak q}=h_{\mathfrak q}(\Delta)-T_{\mathfrak q}
\]

at a prime ideal away from the level.  The following must be proved with one
fixed double-coset normalization:

- the exact function `h_q` on both two-cusp Eisenstein channels;
- annihilation of the residual spectrum as well as the continuous spectrum;
- a rigorous operator-norm bound for `T_q`;
- the functional-calculus Lipschitz constant;
- a certified Hecke defect and non-resonance lower bound for the trial.

No formula has been guessed from the real modular group.  This branch remains
an explicit assumption in the theorem draft.

## Ranked next work

1. **Projected norm certificate** -- cheapest falsification test.  If
   `P_target` retains substantial mass, stay on the symmetry route.  If it
   nearly kills the trial, abandon the shortcut before building the atlas.
2. **Interval first derivatives of overlap defects** -- required by any
   honest smooth-gluing proof.
3. **Explicit finite-type Humbert partition and transition enumeration** --
   converts the abstract lemma into numbers `b0,b1`.
4. **Cancellation-preserving box model** -- essential numerically; refinement
   of the current absolute derivative majorant will not bridge the gap.
5. **General Hecke lemma** -- important for newforms, but not on the critical
   path if the target projector is certified nonzero.

## Sources checked

- J. S. Friedman, *The Selberg trace formula and Selberg zeta-function for
  cofinite Kleinian groups with finite-dimensional unitary representations*,
  Chapter 3 (local copy: `bianchiselberg-refs/friedman_thesis.pdf`).  It
  defines singular subspaces and vector-valued Fourier channels.
- J. S. Friedman, *Analogues of the Artin factorization formula...* for
  induced representations of cofinite Kleinian groups:
  <https://arxiv.org/abs/math/0702030>.
- E. Brenner and F. Spinu, *Artin formalism for Selberg zeta functions of
  co-finite Kleinian groups*:
  <https://www.numdam.org/articles/10.5802/jtnb.657/>.
- A. Booker, A. Strombergsson, and A. Venkatesh, *Effective computation of
  Maass cusp forms*, especially the quasimode/cuspidalizer construction:
  <https://math.stanford.edu/~akshay/research/bsv.pdf>.
- W. Ballmann and P. Polymerakis, *On the differential form spectrum of
  geometrically finite orbifolds*, especially Remark 3.14 on essential
  self-adjointness over complete orbifolds:
  <https://arxiv.org/abs/2011.13304>.
- H. Then, Picard symmetry classes and eigenvalue table (local copy:
  `bianchiselberg-refs/then_thesis.pdf`).

## Artifacts

- `theorem_DK_sixcopy.tex` -- theorem draft
- `six_copy_DK_algebra.py` -- exact permutation/holonomy verifier
- `six_copy_DK_algebra_result.json` -- machine-readable exact result
- `six_copy_target_symmetry.py` -- floating target-class diagnostic
- `six_copy_target_symmetry_result.json` -- diagnostic record
- `track_b_projected_mass_arb.py` -- interval Track-B cusp-witness certificate
- `track_b_projected_mass_arb_result.json` -- canonical 192-bit/512-segment result
