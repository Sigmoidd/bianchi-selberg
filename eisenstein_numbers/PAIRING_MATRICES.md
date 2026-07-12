# Pairing matrices for F₃ / P₃ (g ∈ Γ)

**Purpose.** Freeze explicit elements of Γ = PSL(2, ℤ[ω]) used as
Poincaré side-pairing generators and as residue-gluing labels for the
multi-copy CR assembly. Closes the paper residual “pairing lemma /
computational g∈Γ” for Theorems A–B in `PROOF.md`.

**Reproduce:** `python -u pairing_matrices.py`

---

## 1. Generators as matrices in SL(2, O)

O = ℤ[ω], ω² + ω + 1 = 0. All four matrices have det = 1 ∈ O\* and
descend to PSL(2, O) = Γ.

| Name | Matrix | Möbius on ℂ ∪ {∞} | Order in PSL |
|------|--------|-------------------|--------------|
| T₁ | \(\begin{pmatrix} 1 & 1 \\ 0 & 1 \end{pmatrix}\) | z ↦ z+1 | ∞ (parabolic) |
| T_ω | \(\begin{pmatrix} 1 & \omega \\ 0 & 1 \end{pmatrix}\) | z ↦ z+ω | ∞ (parabolic) |
| U | \(\begin{pmatrix} \omega^2 & 0 \\ 0 & \omega \end{pmatrix}\) | z ↦ ωz | 3 |
| S | \(\begin{pmatrix} 0 & -1 \\ 1 & 0 \end{pmatrix}\) | z ↦ −1/z | 2 |

Code freeze: `geometry_fund.GEN` (complex entries) and
`pairing_matrices.GEN_O` (exact Eisenstein entries).

**Inverses (same matrices up to units in PSL where noted):**

```text
T₁⁻¹ = [[1, −1], [0, 1]]
T_ω⁻¹ = [[1, −ω], [0, 1]]
U⁻¹  = [[ω, 0], [0, ω²]]     # z ↦ ω² z
S⁻¹  = S
```

---

## 2. Faces of the computational cell and reference maps

Planar domain P₃ ([EGM98, §7.3], [DP20, §2.3]):

```text
LEFT  : x = 0,           0 ≤ y ≤ 1/√3
RIGHT : x = 1/2,    |y| ≤ 1/(2√3)
UP    : y = (1−x)/√3,    0 ≤ x ≤ 1/2
LOW   : y = −x/√3,       0 ≤ x ≤ 1/2
FLOOR : |z|² + t² = 1 over P₃
TOP   : t = Y            (free; cusp ODE)
```

| Generator δ | Edges / faces paired on P₃ | Reference face map φ_δ |
|-------------|----------------------------|-------------------------|
| T₁ | RIGHT ↔ LEFT | height-matched opposite sides of the P₃ strip (dictionary) |
| T_ω | LOW ↔ UP | height-matched opposite slanted sides |
| U | vertical ↔ vertical | φ_U(z,t) = (ωz, t) (nearest vertical face) |
| S | FLOOR ↔ FLOOR | full H³ action of S (sphere inversion) |
| — | TOP | free (Lemma D0) |

Production dictionary: `face_pairings_p3.build_pair_maps`.  
Multi-copy rule (same as Gaussian `CONGRUENCE.md` §3):

```text
copy c, δ-source face  ─glue─  copy j = c · δ⁻¹, δ-target face
via reference map φ_δ in the pulled-back cell coordinates.
```

If j = c the identification is internal and is **relaxed** (Neumann),
exactly as at level 1. If j ≠ c the interface is interior to the
Γ₀-domain and must be **conforming**.

---

## 3. Pairing lemma (statement)

**Lemma (pairing matrices).**  
(1) The four matrices of §1 lie in SL(2, O) ⊂ SL(2, ℂ) and represent
elements of Γ.  
(2) **S:** for every floor face f in the certified mesh dictionary,
S(f) lands on the dictionary partner (within mesh tolerance); S² = id.  
(3) **U:** for every vertical dictionary pair labeled U, the image of
the source centre under z ↦ ωz (or z ↦ ω²z) matches the partner centre
up to mesh tolerance.  
(4) **T₁, T_ω:** as lattice translations they generate (with U) the
cusp stabilizer Γ_∞; the P₃ edge dictionary RIGHT↔LEFT / LOW↔UP is the
combinatorial section used for residue gluing of opposite edges of the
EGM strip (area |T| = √3/6). Their inverses act on ideal vertices by
the user cycles in `GENERATOR_CYCLES.md` (Γ_∞-equivalent).  
(5) Residue right action by δ⁻¹ on ℙ¹(O/𝔭) is bijective with
S² = U³ = id (checked for N ∈ {3,7,13}).

Parts (1)–(3) and (5) are machine-checked by `pairing_matrices.py` and
`check_generators.py`. Part (4) records the EGM/Swan section convention
(greedy Γ_∞ reduction is multi-valued; user cycles fix the polyhedron
section).

---

## 4. Ideal-vertex cycles (inverses)

Ground truth (user):

| g⁻¹ | cycles on {∞,0,1,ω,ω²} |
|-----|------------------------|
| T₁⁻¹ | (∞)(0 1)(ω ω²) |
| T_ω⁻¹ | (∞)(0 ω)(1 ω²) |
| U⁻¹ (user label) | (∞)(0)(1 ω ω²) |
| S⁻¹ | (∞ 0)(1)(ω ω²) |

Note: user U⁻¹ cycle equals our **forward** U (z ↦ ωz); matrix U⁻¹ is
the opposite 3-cycle. Same cyclic group ⟨U⟩; residue gluing already
uses δ⁻¹.

---

## 5. What the machine cert uses

| Layer | Object |
|-------|--------|
| Residue | exact perms from matrices mod 𝔭 (`residue_omega.gluing_perms`) |
| Geometry | edge dictionary φ_δ (`face_pairings_p3`) |
| Linear algebra | conforming CR glue when j ≠ c; Neumann when j = c |
| Tops | never in pair maps ⇒ Lemma D0 |

---

## 6. Relation to the Gaussian table

Compare `independent_exclusion/CONGRUENCE.md` §3:

| Gaussian δ | ω analogue |
|------------|------------|
| T₁ = (1,1;0,1) | T₁ same |
| R = rotation by i | U = rotation by ω (order 3) |
| T_i R | lattice conjugates of T₁ / T_ω |
| S = (0,−1;1,0) | S same |

The ω cell P₃ is the unit-folded cusp section (area covol(O)/3), so
opposite-edge dictionary maps for T₁/T_ω are the strip section of the
hexagonal lattice pairings rather than a unit-width parallelogram map.
