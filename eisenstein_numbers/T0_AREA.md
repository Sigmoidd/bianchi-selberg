# Proof: cross-section areas for Γ₀(𝔭) ⊂ PSL(2, ℤ[ω])

**Proposition.** Let O = ℤ[ω], ω = e^{2πi/3}, Γ = PSL(2, O), and let
𝔭 ⊂ O be a **prime ideal** with N(𝔭) = q < ∞ such that O/𝔭 is a field
(equivalently: 𝔭 is ramified or split, so q is prime in ℤ; the first cases
are q ∈ {3,7,13,…}). Write

```text
Γ₀(𝔭) = { [a b ; c d] ∈ Γ  :  c ∈ 𝔭 }.
```

Then Γ₀(𝔭)\ℍ³ has **exactly two cusps** (∞ and 0), and the Euclidean
areas of the orbifold cusp cross-sections (after the standard unit fold)
are

```text
|T_∞| = √3 / 6,
|T_0| = q · √3 / 6 = N(𝔭) · |T_∞|.
```

Consequently the multi-cusp zero-mode weights in the independent-exclusion
criterion are

```text
β_∞(s) = (1−s) / (|T_∞| Y²),
β_0(s)  = (1−s) / (|T_0|  Y²) = (1−s) / (N(𝔭) |T_∞| Y²),
κ_c(s)  = 1 / ((1+s) Y²)     (same Y for both cusps).
```

This is the ℤ[ω]-analogue of `independent_exclusion/CONGRUENCE.md` §2 and
`0cuspchecks.md` C1–C5 (Gaussian O = ℤ[i]).

**Used by:** `cusps_omega.py`, `congruence_omega_proto.py`.

---

## 0. Level-1 data for Γ = PSL(2, O)

Identify ℂ ≅ ℝ². As a lattice,

```text
O = ℤ ⊕ ℤ ω ⊂ ℂ,   covol_ℝ²(O) = |Im(ω)| = √3 / 2.
```

The unit group is O\* = {±1, ±ω, ±ω²} (order 6). In **PSL**, −I is
trivial, so the image of O\* in the cusp stabilizer acts on ℂ by
**rotations of order 3** (multiplication by ω). Write

```text
U = ⟨ z ↦ ω z ⟩  ⊂ Isom(ℂ),   |U| = 3.
```

The parabolic stabilizer is the semidirect product

```text
Γ_∞ ≅ O ⋊ U
```

(upper-triangular matrices with unit diagonal, modulo ±I). A fundamental
domain for O ⋊ U acting on ℂ is the EGM planar domain P_3 of area

```text
|T| := area(P_3) = covol(O) / |U| = (√3/2) / 3 = √3 / 6.
```

(Cf. `GEOMETRY.md`, `geometry_fund.area_P3`, [EGM98, §7.3], [DP20, §2.3].)
This is the level-1 orbifold cross-section: **|T_∞| for Γ equals |T|**.

Shimizu: for γ ∈ Γ with bottom-left entry c ≠ 0, N(c) ≥ 1 (minimal
nonzero ideal norm in O is 1), so {y > Y} with Y ≥ 1 is precisely
invariant under Γ except through Γ_∞ ([DESIGN.md] Lemma 0, same proof
with N(c) ≥ 1 for O = ℤ[ω]).

---

## 1. Cosets and two cusps

Right cosets Γ₀(𝔭)\Γ are labeled by bottom rows of matrix representatives:

```text
[a b; c d]  ↦  (c : d) ∈ ℙ¹(O/𝔭).
```

When O/𝔭 ≅ 𝔽_q is a field, |ℙ¹(𝔽_q)| = q + 1, with transversal

```text
(0 : 1)   and   (1 : x)  for x ∈ 𝔽_q.
```

**Cusp of a coset.** The cusp end of the copy γF sits at γ(∞) = a/c
(in lowest terms). This is Γ₀(𝔭)-equivalent to ∞ if and only if c ∈ 𝔭,
i.e. the label has first coordinate 0 in 𝔽_q. Hence:

| cusp class | labels | count |
|------------|--------|------:|
| ∞ | (0:1) | 1 |
| 0 | (1:x), x∈𝔽_q | q |

So there are **two cusps**, with coset multiplicities (1, q).  
(Verified in code: `cusps_omega.verify_cusp_classes`.)

---

## 2. Cusp ∞: |T_∞| = √3/6

### 2.1 Stabilizer

An upper-triangular element of Γ has c = 0, hence lies in Γ₀(𝔭) for
every 𝔭. Thus

```text
Γ_∞ ∩ Γ₀(𝔭) = Γ_∞ ≅ O ⋊ U.
```

The cusp-∞ cross-section is therefore **identical** to the level-1
orbifold section:

```text
|T_∞| = area(P_3) = √3 / 6.
```

### 2.2 Collar

Truncate the level-1 fundamental domain F at height Y. The cusp-∞ collar
of Γ₀(𝔭) is still P_3 × (Y,∞): the horoball {y > Y} meets its
Γ₀(𝔭)-translates only through Γ_∞ (since Γ₀(𝔭) ⊂ Γ and Lemma 0 applies
to Γ). One coset copy of class ∞ contributes top area |T_∞| to the
trace functional t_∞.

---

## 3. Cusp 0: lattice, fold, and |T_0|

Let σ₀ = S = \(\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\) ∈ Γ, so
σ₀(∞) = 0. Work in the conjugate chart that sends 0 to ∞.

### 3.1 Conjugated stabilizer (lattice + unit fold)

A matrix γ = \(\begin{pmatrix} a & 0 \\ c & d \end{pmatrix}\) ∈ Γ₀(𝔭)
fixes 0 (so b = 0), with ad − bc = ad a unit (here = 1 after scaling into
SL), and **c ∈ 𝔭**. Conjugation by σ₀ (direct matrix multiply, same
pattern as `0cuspchecks.md` C1):

```text
σ₀⁻¹ (a  0 ; c  d) σ₀  =  (d  −c ; 0  a).
```

(Up to the central ±I of PSL.) This is upper triangular with:

- **translations** z ↦ z + τ where τ runs through 𝔭 (as c runs through 𝔭,
  absorbing units into the diagonal), hence

```text
Λ₀ = 𝔭  ⊂ ℂ   (ideal lattice).
```

  As a full-rank lattice,  
  `covol_ℝ²(𝔭) = N(𝔭) · covol(O) = q · √3 / 2`.

- **unit diagonals** (d, a) with ad = 1. Taking a = ω (so d = ω², since
  ω·ω² = 1) gives multiplication by ω² on ℂ after conjugation — an
  element of order 3 in PSL. Thus the **same order-3 fold U** is present
  at cusp 0 for every Γ₀(𝔭).

**Cross-section area:**

```text
|T_0| = covol(Λ₀) / |U| = (q · √3 / 2) / 3 = q · √3 / 6 = N(𝔭) · |T_∞|.
```

This is the claimed formula. ∎

### 3.2 Fold-compatible tiling of T_0

Let T ⊂ ℂ be a fundamental domain for O ⋊ U of area |T_∞| (e.g. P_3).
Let R be a complete set of residues for O/𝔭 (so |R| = q). Claim:

```text
T_0 :=  ⊔_{r ∈ R}  (T + r)
```

is a fundamental domain for Λ₀ ⋊ U = 𝔭 ⋊ U.

**Proof** (semidirect product, as in `0cuspchecks.md` C2). An element of
O ⋊ U acts by z ↦ u z + λ with u ∈ U, λ ∈ O. It lies in
(𝔭 ⋊ U)·τ_r (τ_r = translation by r) iff

```text
λ − u r ∈ 𝔭  ⟺  r ≡ u⁻¹ λ   (mod 𝔭).
```

As r runs through a complete set of residues of O/𝔭 and u is a unit
(automorphism of O/𝔭), there is **exactly one** such r for each
(λ, u). The fold U preserves 𝔭 (ideal) and acts on O/𝔭, so the
decomposition is U-equivariant. Areas add:

```text
area(T_0) = q · area(T) = q · |T_∞| = |T_0|,
```

matching §3.1. ∎

### 3.3 Chimneys at common height Y (Fact B)

For each residue class x ∈ 𝔽_q choose a lift \(\tilde x ∈ O\) and set

```text
γ_x = \begin{pmatrix} 0 & -1 \\ 1 & \tilde x \end{pmatrix}
```

(bottom row (1 : x)). Together with the identity (bottom row (0:1)), these
exhaust ℙ¹(𝔽_q). In PSL one has S⁻¹ = S, and

```text
S · γ_x  ≡  \begin{pmatrix} 1 & \tilde x \\ 0 & 1 \end{pmatrix}
```

(a **height-preserving** translation). Hence in cusp-0 coordinates the
class-0 chimneys are

```text
(T + \tilde x) × (Y, ∞),  x ∈ 𝔽_q,
```

which assemble exactly to T_0 × (Y, ∞) at the **same** truncation height Y
as cusp ∞. ∎

### 3.4 Trace assembly identity

Each of the q coset copies of class 0 has a reference top face of Euclidean
area |T_∞|. The change of coordinates to the cusp-0 collar is a
height-preserving isometry (translation), so

```text
t_0(1) := Σ_{class-0 copies} ∫_{top} 1 dx  =  q · |T_∞|  =  |T_0|.
```

Similarly t_∞(1) = |T_∞|. This is the check in
`congruence_omega_proto.py` (N=3: t_∞=√3/6, t_0=√3/2).

---

## 4. Zero-mode defect and β_α

On each cusp collar T_α × (Y,∞), the zero Fourier mode of an
L² eigenfunction with λ = 1−s² contributes
(cf. DESIGN.md Lemma 2, |T|-parametrized)

```text
− |T_α| (1−s) c_α² / Y²,   c_α = mean of u(·,Y) on T_α.
```

With t_α(v) := ∫_{T_α} v dx = |T_α| c_α,

```text
defect_α = − (1−s) t_α(v)² / (|T_α| Y²).
```

Hence the pencil coefficients are exactly

```text
β_α(s) = (1−s) / (|T_α| Y²),
```

i.e.

```text
β_∞(s) = (1−s) / (|T_∞| Y²),
β_0(s)  = (1−s) / (N(𝔭) |T_∞| Y²).
```

Orthogonality to constants gives a single rank-one constraint with

```text
κ_c(s) = 1/((1+s) Y²)
```

on t_∞ + t_0 (same Y on both collars). ∎

---

## 5. Nonzero-mode bound at cusp 0 (float / design check)

The dual lattice of Λ₀ = 𝔭 scales as N(𝔭)^{−1/2} relative to the dual of
O. With level-1 dual length |μ|_∞ ≳ 2/√3 (hexagonal),

```text
|μ|_0  ≳  (2/√3) / √q,
4π² Y² / |μ|_0²  ≳  4π² Y² · 3 / (4 q) = 3 π² Y² / q.
```

At Y = 1.25 and q ≤ 13 this is ≫ 1, so Lemma 1 (nonzero modes ≥ 0)
applies on both cusps for the first ladder rungs. (Folding to the
U-invariant subspace only shrinks the test space further.)

---

## 6. Volume sanity check

Fundamental domain volume multiplies by the index:

```text
vol(Γ₀(𝔭)\ℍ³) = (q+1) · vol(Γ\ℍ³).
```

Truncation:

```text
vol(D_Y) = (q+1) vol(K_Y)
         + |T_∞|/(2Y²) + |T_0|/(2Y²)
         = (q+1)( vol(F) − |T|/(2Y²) )
         + (1+q) |T| / (2Y²)
         = (q+1) vol(F),
```

using |T_∞| = |T|, |T_0| = q|T|, vol(K_Y) = vol(F) − |T|/(2Y²). Checks. ∎

---

## 7. Scope and non-claims

| Claim | Status |
|-------|--------|
| \|T_∞\| = √3/6 for Γ and for Γ₀(𝔭) | **proved** (§0–2) |
| \|T_0\| = N(𝔭)·√3/6 for prime 𝔭 with O/𝔭 a field | **proved** (§3) |
| Two cusps with coset counts (1,q) | **proved** (§1) |
| β_α, κ_c as above | **proved** (§4) |
| Inert primes (O/𝔭 ≅ 𝔽_{p²}, e.g. p=2) | **not** covered — residue is not a prime field; redo lattice/fold |
| Composite square-free 𝔫 via CRT | **not** covered — product of local patterns (`P1_ACTION.md`) |
| Interval certification of Γ₀ | independent of this note; needs Arb/Rump on glued FE space |

---

## 8. References

1. Elstrodt–Grunewald–Mennicke, *Groups Acting on Hyperbolic Space*, §7.2–7.3 (F_d, cusp structure).
2. Dória–Paula, arXiv:1910.03148, §2.3 (explicit P_3).
3. In-repo Gaussian parallel: `independent_exclusion/CONGRUENCE.md` §1–2, `0cuspchecks.md` C1–C5.
4. Level-1 reduction: `independent_exclusion/DESIGN.md` Lemmas 0–2 (|T|-form).
5. Code: `cusps_omega.py`, `congruence_omega_proto.py`, `residue_omega.py`.
