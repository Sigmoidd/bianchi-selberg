# Stage 8 — Congruence ladder (PSL(2, ℤ[ω]))

**Status:** foundation started (`residue_omega.py`). No eigenvalue certificate yet.  
Parallel to `independent_exclusion/CONGRUENCE.md` for Γ₀(𝔫) ⊂ PSL(2, ℤ[i]).

## Target

> For prime ideals 𝔭 ⊂ ℤ[ω] (then square-free products),  
> Γ₀(𝔭)\ℍ³ has **no Laplace eigenvalue in (0,1)**.

Method: Reference-Cell Principle + multi-cusp G1 + Arb/Rump on glued CR space  
(same as Thms 1–4 for ℤ[i], different residue geometry).

## Prerequisites

| # | Item | Status |
|---|------|--------|
| 1 | Level-1 FEM cert PSL(2,ℤ[ω]) | **done** — `PROOF.md`, `cert_omega.py` |
| 2 | Ideal / residue arithmetic | **done** — `residue_omega.py` (N=3,7,13) |
| 3 | ℙ¹(R) right action + gluing perms | **done** — field case |
| 4 | Cusp classification Γ₀(𝔭) | **table** — `cusps_omega.py` (prove \|T_0\| in paper) |
| 5 | Multi-cusp criterion (β_α, \|T_α\|) | **coded** — same shape as Gaussian §2 |
| 6 | Float CR on N(𝔭)+1 copies | **float PASS** N=3 @ mesh 6×3 |
| 7 | Interval cert (m3p-style) | **CERTIFIED N=3** — `cert_omega_p.py` 8/8 |

## First rungs (proposed order)

| Rung | 𝔭 | N(𝔭) | index \|ℙ¹\| | Residue | Notes |
|------|---|------|-------------|---------|--------|
| 0 | level 1 | 1 | 1 | — | **certified** (`cert_omega.py`) |
| 1 | (1−ω) | 3 | 4 | 𝔽₃ | ramified; smallest index |
| 2 | π \| 7 | 7 | 8 | 𝔽₇ | split p≡1 mod 3 |
| 3 | π \| 13 | 13 | 14 | 𝔽₁₃ | parallel to Gaussian N=13 |
| later | (2) inert | 4 | 5 | 𝔽₄ | needs quadratic residue ring |
| later | products | CRT | … | ∏ 𝔽_{qᵢ} | P1_ACTION pattern |

Norm formula: N(a+bω)=a²−ab+b².  
Units: 6. Index for prime 𝔭 with field O/𝔭: **N(𝔭)+1**.

Reproduce inventory:

```powershell
python -u residue_omega.py
```

## What carries over from ℤ[i]

| Ingredient | Carry? |
|------------|--------|
| Multi-cusp Lax–Phillips criterion | Yes (per-cusp \|T_α\|, β_α, κ_c) |
| Theorem G1 / G1𝔭 structure | Yes (level-1 cell = EGM K_Y) |
| ℙ¹(R) right action | Yes — field case in `residue_omega.py` |
| Rump equilibration | Yes (`m3p_certify.rump_certify_inplace`) |
| Reference-Cell Principle | Yes — one P₃ cell, copies = index |
| Face pairings T1,R,TiR,S (Gaussian) | **No** — use T1, Tw, U, S (hex) |
| Residue rings with i↦… | **No** — ω↦w_img in 𝔽_q |

## Gluing generators (level-1 faces)

From `geometry_fund.GEN` / `GEOMETRY.md` §4:

| δ | matrix | role |
|---|--------|------|
| T₁ | [[1,1],[0,1]] | vertical wall u-pair |
| T_ω | [[1,ω],[0,1]] | lattice translation |
| U | [[ω²,0],[0,ω]] | order-3 unit (z↦ωz) |
| S | [[0,−1],[1,0]] | floor |

Right action on ℙ¹(R): copy p glues via δ to copy p·δ⁻¹ (same convention as
`CONGRUENCE.md` §3). Internal fixed points of the perm → Neumann (relaxed).

## Multi-cusp criterion (proved for field primes)

**Proof:** [`T0_AREA.md`](T0_AREA.md) (ℤ[ω] parallel of `0cuspchecks.md` C1–C5).

For prime 𝔭 with N(𝔭)=q and O/𝔭 a field, ℙ¹ partitions as **(1 : q)**
for cusps (∞ : 0) — verified in `cusps_omega.verify_cusp_classes`.

| cusp | cosets | \|T_α\| | β_α |
|------|--------|---------|-----|
| ∞ | 1 | √3/6 | (1−s)/(\|T_∞\| Y²) |
| 0 | q | **q·√3/6** | (1−s)/(\|T_0\| Y²) |

κ_c = 1/((1+s)Y²) on t_∞+t_0 (rank-1 constraint).

Outline: Λ₀ = 𝔭 has covol q·√3/2; order-3 unit fold still present after
conjugating by S ⇒ |T_0| = covol(Λ₀)/3 = q·√3/6. Chimneys T+r tile T_0
at common height Y.

### Float result (rung 1)

```text
python -u congruence_omega_proto.py 3 6 3
# N=3, index 4, connected, t_∞/t_0 exact, min μ ≈ 0.90 > 0
```

Face pair maps: **EGM P_3 edge dictionary** (`face_pairings_p3.py`) —

| generator | edges paired |
|-----------|----------------|
| T1 | RIGHT ↔ LEFT |
| Tw | LOW ↔ UP |
| U | vertical faces via z↦ωz |
| S | FLOOR ↔ FLOOR (sphere inversion) |

Float @ 6×3: min μ ≈ **1.12** (improved vs earlier tags).

## Implementation plan

```text
residue_omega.py           # ring, P1, gluing perms           [done]
cusps_omega.py             # cusp areas / β / modes           [done]
T0_AREA.md                 # proof |T_0|=N|T_∞|               [done]
face_pairings_p3.py        # EGM edge pairing dictionary      [done]
congruence_omega_proto.py  # float multi-copy CR + pencil     [float PASS N=3]
CERT_OMEGA_P.md            # cert design / G1p / memory       [done]
cert_omega_p.py            # Arb + Rump                       [**CERTIFIED N=3**]
```

### Certificate freeze (N=3, 2026-07-12)

```text
mesh P3 6×3, Y=1.25, ν*=1.05
(θ, θ₂, α, θ₄, ρ̃) = (0.5, 0.9, 0.2, 0.5, 9.0)
κ = Lemma I1, Taylor p=5, arb prec 128, Rump SAS+per-row
reproduce: python -u cert_omega_p.py 3 6 3
```

## Non-goals (until next rungs)

- N=7,13 certificates; CRT composites; inert 𝔽₄.
- Pairing matrices frozen: `PAIRING_MATRICES.md` / `pairing_matrices.py` PASS.
