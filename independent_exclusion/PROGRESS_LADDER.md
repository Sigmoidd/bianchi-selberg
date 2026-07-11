# Ladder progress log (worktree `np9` only)

All code and runs live under
`C:\Users\Admin\.grok\worktrees\bearings-bianchi-selberg\np9\independent_exclusion`.
The main tree at `Documents\Bearings\bianchi-selberg` was not modified.

## Goal

Certify: no Laplace eigenvalues in (0,1) on Γ₀(𝔭)\ℍ³ for
- 𝔭 = (3) (inert, N𝔭 = 9, index 10)
- 𝔭 = (3+2i) (split, N𝔭 = 13, index 14)

Plan of record: `HANDOFF_LADDER.md` / `CONGRUENCE.md` §8.

---

## Done: step 1 — residue ring + M𝔭0 margin (N𝔭 = 9)

### Code changes

**`congruence_prototype.py`**
- Replaced hardcoded 𝔽₅ arithmetic with `GaussianResidueRing` + `set_level(level)`.
- Supported levels:
  | level | ring | N(p) | i-image |
  |-------|------|------|---------|
  | `(2+i)` | 𝔽₅ prime field | 5 | 3 |
  | `(3)` | 𝔽₉ as pairs (a,b) mod 3, i²=−1 | 9 | (0,1) |
  | `(3+2i)` | 𝔽₁₃ prime field | 13 | 8 |
- Generators: T1, R_rot=(i,0;0,−i), T_iR=(1,i;0,1)·R_rot, S over R.
- Asserts: permutations, round-trip act, involutions M²≡±I in PSL for R/TiR/S, cusp classes {1, N(p)}, connectivity.
- Assembly: `NC = NP + 1`; checks `t_∞(1)=1/2`, `t₀(1)=N(p)/2`, connected dof graph.
- CLI: `python congruence_prototype.py gluing` |
  `python congruence_prototype.py margin "(3)"` |
  `python congruence_prototype.py validate "(3)"`.

**`m3p_certify.py`**
- Live NP/LEVEL via `congruence_prototype` module (after `set_level`).
- Lemma D0 constant: `σ_h = √(N(p)+1) · α_h^ref` (was √6).
- β₀ uses live NP; CLI: `python m3p_certify.py coarse "(3)"`.

No changes to: reference mesh builder, arb element enclosures, Rump core, window formulas (other than NP / √(NP+1)).

### Gluing assert suite (all levels)

```
gluing OK  level=(2+i)  N(p)=5   copies=6   cusps inf=1 zero=5
gluing OK  level=(3)    N(p)=9   copies=10  cusps inf=1 zero=9
gluing OK  level=(3+2i) N(p)=13  copies=14  cusps inf=1 zero=13
```

F₅ generator matrices bit-identical to the previous hardcode
(`R=[[3,0],[0,2]]`, `TiR=[[3,1],[0,2]]`, `S=[[0,4],[1,0]]`).

### M𝔭0 float margin — Γ₀(3), reference 6×3×3

Command:
```text
python -u congruence_prototype.py margin "(3)"
```

| quantity | value |
|----------|-------|
| copies | 10 |
| tets/copy | 1296 |
| global CR dofs | 26 532 |
| 1′M1 | 1.382075 (≈ 10 · 0.1382) |
| t_∞(1) | 0.5 |
| t₀(1) | 4.5 (= N(p)/2) |
| \|Q·1\| | ~4e-16 |
| connected components | 1 |
| **min μ over λ∈(0,1)** | **≈ 2.399** |

Coarse 4×2×2 cross-check: min μ ≈ 2.45 (same order).

Handoff expected μ ≈ 3–3.5; measured **2.40** is lower than the
sublinear extrapolation (7.3 → 4.4 at N=5) but **well above the ~1.5
audit floor**. Assembly invariants match §4 of the handoff → proceed
to certification.

Reproduce:
```powershell
cd C:\Users\Admin\.grok\worktrees\bearings-bianchi-selberg\np9\independent_exclusion
python -u congruence_prototype.py gluing
python -u congruence_prototype.py margin "(3)"
```

---

## Done: step 2 — certify N𝔭 = 9 (Theorem 3)

**Result: CERTIFIED** — no eigenvalues in (0,1) for Γ₀(3)\ℍ³.
Wall time ≈ 38 min (2026-07-11). Written as **Theorem 3** in `PROOF.md`.

### Parameter search note

Default Np=5-style grid preferred high ρ̃ and failed the pencil filter
(best min-eig ≈ 0.90). Expanded grid found lower ρ̃ feasible; frozen:

| param | value |
|-------|-------|
| (θ, θ₂, α, θ₄, ρ̃) | **(0.5, 0.9, 0.15, 0.8, 1.5)** |
| ν* | 1.02 |
| mesh | 6×3×3, 26 532 dofs |
| float pencil (search windows) | ≈ 1.13 |

### Per-window Rump certificate

| window s | c_S>0 | c_e≥d_e | PSD | c_e | d_e |
|----------|-------|---------|-----|-----|-----|
| [0.000,0.125] | True | True | True | 0.836 | 0.222 |
| [0.125,0.250] | True | True | True | 0.837 | 0.219 |
| [0.250,0.375] | True | True | True | 0.840 | 0.210 |
| [0.375,0.500] | True | True | True | 0.845 | 0.196 |
| [0.500,0.625] | True | True | True | 0.851 | 0.175 |
| [0.625,0.750] | True | True | True | 0.859 | 0.149 |
| [0.750,0.875] | True | True | True | 0.869 | 0.117 |
| [0.875,1.000] | True | True | True | 0.881 | 0.079 |

Reproduce:
```powershell
cd C:\Users\Admin\.grok\worktrees\bearings-bianchi-selberg\np9\independent_exclusion
python -u -c "from m3p_certify import certify; print(certify(6,3,3,params=(0.5,0.9,0.15,0.8,1.5),level='(3)'))"
```

---

## Step 3 — N𝔭 = 13: CERTIFIED (Theorem 4)

**Result: CERTIFIED** after power-of-two diagonal equilibration + per-row
radius Rump (2026-07-11). Wall time ≈ 93 min. See `PROOF.md` Theorem 4.

### Key fix (why unscaled Rump failed)

Conservative midpoint A was SPD with λ_min ≈ 3·10⁻⁵, but a **uniform**
max-row-sum radius shift ≈ 1.5·10⁻⁴ exceeded that margin. Equilibration
s_i = 2^round(−½ log₂ a_ii) plus **per-row** r_i on the diagonal clears
all 16 windows. Float pencil ~1.13 was never the blocker.

### Frozen certificate (Attempt B + scaling)

| field | value |
|-------|--------|
| mesh | 8×3×2, 33 176 dofs |
| ν* / NWIN | 1.001 / 16 |
| params | (0.3, 0.98, 0.1, 0.85, 0.5) |
| all 16 windows | PSD **True** |
| c_e range | 0.870 … 0.907 |
| d_e range | 0.074 … 0.485 |

Reproduce:
```powershell
python -u -c "import m3p_certify as m; m.NU_STAR=1.001; m.NWIN=16; print(m.certify(8,3,2,params=(0.3,0.98,0.1,0.85,0.5),level='(3+2i)'))"
```

---

## Appendix — earlier N𝔭=13 diagnostics (superseded by Theorem 4)

### Float margin M𝔭0 (assembly OK everywhere)

| mesh | dofs | min μ | notes |
|------|------|-------|-------|
| 6×3×2 | 24 888 | **≈ 1.445** | preferred handoff mesh |
| 6×3×3 | 37 044 | ≈ 1.454 | same order |
| 8×4×2 | 44 192 | ≈ 1.436 | μ not helped by refine |
| 8×4×3 | 65 776 | ≈ 1.446 | μ mesh-converged ~1.45 |

Checks: 14 copies, t_∞(1)=0.5, t₀(1)=6.5, connected. Handoff expected
μ ≈ 2.5–3; measured **~1.45** (at the ~1.5 audit floor). Steeper decay
than 7.3 → 4.4 → 2.4 suggested.

### Certify attempts (all np9; main untouched)

| attempt | mesh | ν* | NWIN | params (θ,θ₂,α,θ₄,ρ̃) | float pencil | Rump |
|---------|------|-----|------|------------------------|--------------|------|
| A | 6×3×2 | 1.005 | 16 | (0.35, 0.95, 0.1, 0.85, 0.5) | ≈ 1.04 | **PSD False** w0 (c_e=0.85>d_e=0.53) |
| B | 8×3×2 | 1.001 | 16 | (0.3, 0.98, 0.1, 0.85, 0.5) | ≈ 1.13 | **PSD False** w0–w1 (c_e≈0.87>d_e≈0.48) |

Attempt B used planar refinement to cut S_Q (0.061→0.044) and recover
the same float-pencil regime that worked for Theorem 3; Rump still
rejects the mid/rad matrix (Cholesky of diagonally reduced A fails).
Run B stopped after 2 windows to avoid ~2 h of doomed Choleskys
(dense A ≈ 8.8 GB, ~10 min/window).

### Bound analysis (why this is a ceiling)

Two different bounds:

1. **Theory (mode / collar):** N𝔭 < 4π²Y² ≈ 61.7 — N𝔭=13 is fine.
2. **Certificate headroom (what we hit):**  
   `μ_h ≈ 1.45` after G1𝔭 inflation (ν*, c_Q&lt;1 from S_Q, ρ̃, window
   dk, Rump ε) must stay &gt; 1 for PSD. Empirically the inflated form
   sits at float pencil ~1.04–1.13; Rump needs more slack than LOBPCG
   reports. **This is a method bound, not evidence of eigenvalues.**

Margin ladder: 7.3 (lv1) → 4.4 (N=5) → 2.4 (N=9, certified) → **1.45
(N=13, not certified)**. N=9 was the last comfortable rung.

### What would be needed for Theorem 4 (out of current budget)

- More planar resolution with memory for Rump (8×4×3 ~66k dofs ≈ 35 GB
  dense — beyond 16 GB), or out-of-core blocked Cholesky (handoff
  fallback (c), hours), and/or
- Sharper theory constants (S_Q / window lemmas) — HANDOFF §7 out of
  scope, and/or
- Slightly taller Y (changes collar setup; not free).

**Superseded:** N𝔭=13 is certified (Theorem 4) once Rump uses equilibration
+ per-row radii. The μ≈1.45 “audit floor” was not a hard certificate ceiling.

---

## Ladder closeout (np9)

| Step | Status |
|------|--------|
| 1. Residue ring + M𝔭0 | done |
| 2. Certify N𝔭=9 | **Theorem 3** |
| 3. Certify N𝔭=13 | **Theorem 4** |
| 4. PROOF + status | done |
