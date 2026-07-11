# Right action on ℙ¹(ℤ[i]/𝔫)

Companion to `CONGRUENCE.md` §3 and `congruence_prototype.py`.
Describes the combinatorial core for **general** principal level 𝔫,
beyond the prime ladder (Theorems 2–4).

## 1. Setup

Let 𝔫 ⊂ ℤ[i] be a nonzero ideal. Write R = ℤ[i]/𝔫 (finite commutative
ring). Cosets Γ₀(𝔫)\PSL(2,ℤ[i]) are labeled by bottom rows of matrix
representatives:

```text
γ = ( a  b ; c  d )  ↦  (c : d) ∈ ℙ¹(R)
```

**Unimodular pairs** U(R) = { (c,d) ∈ R² : cR + dR = R }.

**Projective line** ℙ¹(R) = U(R) / R^×, where λ·(c,d) = (λc, λd).

When R is a field, |ℙ¹(R)| = |R| + 1 and a complete set of reps is

```text
(0:1)  and  (1:d) for d ∈ R.
```

## 2. Right action (formula unchanged)

For M = ( a b ; c' d' ) ∈ SL(2,R) (PSL: identify ±M),

```text
(c:d) · M  :=  class of  (c,d) M  = (c a + d c',  c b + d d').
```

This is the action used by gluing: copy p's δ-source face meets copy
(p · δ⁻¹)'s δ-target face (`CONGRUENCE.md` §3; `0cuspchecks.md` C7).

Properties:

- Well-defined on R^×-orbits.
- SL(2,R) preserves U(R).
- Face generators T1, R_rot, T_iR, S act by permutations of any
  transversal of ℙ¹(R).

## 3. Field case (prime 𝔭) — already in code

R = ℤ[i]/𝔭 field (or 𝔽₉ for 𝔭=(3)). Normalization:

```text
c ≠ 0  ⇒  (1 : d c⁻¹),   else  (0:1).
```

Cusp class: ∞ iff c ≡ 0, else the unique finite cusp; sizes {1, N𝔭}.

## 4. Square-free composite — CRT product

If 𝔫 = 𝔭₁⋯𝔭_k distinct primes of ℤ[i], then

```text
R ≅ ∏_{j=1}^k R_j,   R_j = ℤ[i]/𝔭_j,
ℙ¹(R) ≅ ∏_{j=1}^k ℙ¹(R_j).
```

Index formula (square-free 𝔫):

```text
[PSL(2,ℤ[i]) : Γ₀(𝔫)] = N(𝔫) ∏_{𝔭|𝔫} (1 + 1/N𝔭) = ∏_j (N𝔭_j + 1).
```

Example: 𝔫 = (5) = (2+i)(2−i), N=25 → index 6·6 = **36**.

**Local embeddings of i** (must satisfy i² = −1 in each R_j):

| prime | N | i image |
|-------|---|---------|
| (2+i) | 5 | 3 |
| (2−i) | 5 | 2 |
| (3)   | 9 | pairs (0,1) in 𝔽₉ |
| (3+2i)| 13| 8 |

**Normalization:** componentwise field normalization; a point is a
k-tuple of local normalized pairs.

**Action:** embed M entrywise into each R_j; act componentwise.

## 5. Cusp labels (combinatorial)

True geometric cusps of Γ₀(𝔫)\ℍ³ are fewer than the naive product of
local {∞, finite} patterns once Γ-identification is taken into account.
For gluing **permutations**, only ℙ¹(R) matters.

Code currently records a **local infinity pattern**

```text
key(p) = (1_{c_j=0})_j
```

for diagnostics. Sizes for 𝔫=(5): (∞,∞):1, (∞,fin):5, (fin,∞):5,
(fin,fin):25. Geometric multi-cusp widths |T_α| for the FEM criterion
are **not** claimed here — that is Facts A–B for composite 𝔫
(`CONGRUENCE.md` §8, out of certificate scope until written).

## 6. What this does *not* include

- Multi-cusp positivity criterion / Theorem G1𝔫
- Rump certificate at composite level
- Non-square-free 𝔫 (p^k) without further local ℙ¹ over R/p^k

## 7. Code map

| Symbol | Implementation |
|--------|----------------|
| Field R | `FieldResidueRing` |
| Product R | `ProductResidueRing` / CRT level `(5)` |
| ℙ¹ + act | methods `p1_points`, `p1_normalize`, `act` |
| Gluing | `build_gluing(level)` |
| Levels | `set_level('(2+i)'|'(3)'|'(3+2i)'|'(5)')` |

Prime regression: gluing asserts for Theorems 2–4 must stay green.
