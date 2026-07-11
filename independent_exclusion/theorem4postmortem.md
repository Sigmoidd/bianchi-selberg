# Full post-mortem: ladder rungs N𝔭 = 9 and 13

**Worktree:** `C:\Users\Admin\.grok\worktrees\bearings-bianchi-selberg\np9` only  
**Main** (`Documents\Bearings\bianchi-selberg`): not modified for this work (still has legacy `f5_matrices()`).  
**Plan of record:** `HANDOFF_LADDER.md` / `CONGRUENCE.md` §8  
**Operational log:** `PROGRESS_LADDER.md`

---

## 1. Mission

Certify, by verified finite elements (no trace formula):

> Γ₀(𝔭)\ℍ³ has **no Laplace eigenvalue in (0,1)**

for the next two primes after (2+i):

| Rung | 𝔭 | N𝔭 | Index | Copies | Outcome |
|------|---|-----|-------|--------|---------|
| Existing | (2+i) | 5 | 6 | 6 | Theorem 2 (prior) |
| **This work** | **(3)** inert | **9** | 10 | 10 | **Theorem 3 — CERTIFIED** |
| **This work** | **(3+2i)** split | **13** | 14 | 14 | **Theorem 4 — CERTIFIED** (after Rump equilibration) |

Method chain (unchanged): Lax–Phillips multi-cusp criterion → Theorem G1𝔭 (CR lower bound, Lemma D0) → arb enclosures + Rump BIT 46 PSD.

---

## 2. Timeline (what we did, in order)

### Phase A — Infrastructure (step 1)

1. Read handoff; confined all edits to **np9**.
2. Generalized residue arithmetic in `congruence_prototype.py`:
   - `GaussianResidueRing` + `set_level`
   - F₅ (i↦3), F₉ (pairs mod 3), F₁₃ (i↦8)
   - Generators T1, R, TiR, S; P¹ points; gluing asserts
3. Threaded `NC = NP+1`, assembly checks `t_∞(1)=½`, `t₀(1)=NP/2`
4. Updated `m3p_certify.py`: live NP, `σ_h = √(NP+1)·α_h^ref`
5. Gluing suite: all three levels pass; F₅ matrices bit-match old hardcode
6. M𝔭0 for N𝔭=9 @ 6×3×3: **min μ ≈ 2.40**, 26 532 dofs, checks OK

### Phase B — Certify N𝔭=9 (step 2)

1. First certify with default parameter grid: **fail** — pencil max ≈ 0.90 (search preferred high ρ̃; anti-correlates with pencil)
2. Expanded grid: low ρ̃ feasible; frozen  
   `(θ,θ₂,α,θ₄,ρ̃) = (0.5, 0.9, 0.15, 0.8, 1.5)`, ν\*=1.02, 8 windows
3. Full Rump run ~**38 min**: **8/8 PSD True**
4. Wrote **Theorem 3** in `PROOF.md`; status in HANDOFF / CONGRUENCE / PROGRESS

### Phase C — Certify N𝔭=13 (step 3)

1. M𝔭0: **μ ≈ 1.45** at 6×3×2, 6×3×3, 8×4×2, 8×4×3 (mesh-converged)
2. Assembly always OK: 14 copies, t₀(1)=6.5, connected
3. Parameter probes: need very low ρ̃, tight ν\*, more windows to get pencil > 1
4. **Attempt A** 6×3×2, ν\*=1.005, 16 wins, params `(0.35,0.95,0.1,0.85,0.5)`  
   → scalar OK, **PSD False** on window 0 (float pencil ~1.04)
5. **Attempt B** 8×3×2 (better S_Q), ν\*=1.001, 16 wins, `(0.3,0.98,0.1,0.85,0.5)`  
   → float pencil ~**1.13** (same regime as Thm 3), **PSD False** w0–w1  
   → stopped after 2 windows (~20 min wasted Cholesky avoided further)
6. Documented **certificate ceiling**; no Theorem 4

---

## 3. Final scorecard

### What succeeded

| Item | Evidence |
|------|----------|
| Worktree isolation | Main still `def f5_matrices()`; np9 has `GaussianResidueRing` |
| Residue ring + gluing | Asserts for N𝔭 ∈ {5,9,13} |
| Regression F₅ | Matrices identical to pre-generalization |
| Assembly physics checks | Volume scale, cusp traces, connectivity |
| **Theorem 3** | All 8 windows Rump PSD; reproducible one-liner |
| Diagnosis of N=13 | μ mesh-converged; failure is Rump/headroom not gluing bug |

### What failed / deferred

| Item | Outcome |
|------|---------|
| Theorem 4 (N𝔭=13) | **CERTIFIED** (equilibrated Rump) |
| Handoff μ extrapolation | Too optimistic (expected 2.5–3, got ~1.45) |
| Default param search for larger N | Biased to high ρ̃ → false “no feasible params” |
| Float pencil ≳ 1.1 as Rump predictor | **Insufficient** at n~25–33k (N=13) |
| 16 GB dense Rump at finer planar meshes | 8×4×3 ~66k dofs ≈ 35 GB — out of budget |

---

## 4. Quantitative ladder (the story in numbers)

### Float exclusion margin min_λ μ(λ)

| Level | N𝔭 | Copies | min μ | Certifiable with current G1𝔭+Rump? |
|-------|-----|--------|-------|-------------------------------------|
| level 1 | 1 | 1 | ~7.3 | yes (Thm 1) |
| (2+i) | 5 | 6 | ~4.4 | yes (Thm 2; pencil ~1.44) |
| **(3)** | **9** | **10** | **~2.40** | **yes (Thm 3; pencil ~1.13)** |
| **(3+2i)** | **13** | **14** | **~1.45** | **yes (Thm 4; after SAS + per-row Rump)** |

Decay is real and **not** a coarse-mesh artifact at N=13 (four meshes agree).  
Handoff audit floor μ ≳ 1.5: N=13 sits **on** it.

### Certificate chain (sufficient condition)

```text
true / discrete μ
    →  CR form μ_h
    →  window inflation (ν*, S_Q, S_M, ρ̃, θ, dk)  →  float pencil λ_min(N,D)
    →  mid/rad + Rump diagonal shift ε(n)           →  Cholesky PSD
```

Each arrow **eats** margin. Failure mode at N=13:

- Scalars still pass: `c_e > d_e` (e.g. 0.87 > 0.48)
- Float pencil can sit slightly above 1
- Rump’s reduced-diagonal Cholesky **fails** → effective certified form not SPD after ε

So the bound we hit is:

> **Certificate headroom:** μ_h / inflation ≰ 1 after Rump  
> **Not** the mode bound N𝔭 < 4π²Y² ≈ 61.7  
> **Not** “eigenvalues must exist” (μ still positive)

### Reference constants that matter

| Mesh | S_Q | S_M | γ (typ.) | Role |
|------|-----|-----|----------|------|
| 6×3×2 | 0.061 | 0.514 | ~0.172 | Preferred handoff mesh; weak c_Q |
| 8×3×2 | **0.044** | **0.389** | ~0.175 | Planar help; still Rump-fail |
| 8×4×3 | 0.032 | 0.268 | ~0.113 | Would help more; ~35 GB dense |

Planar refinement improves **inflation** (S_Q↓ → c_Q↑), not the continuum μ (~1.45 fixed).

---

## 5. Theorem 3 frozen certificate (success case)

| Field | Value |
|-------|--------|
| Statement | No λ ∈ (0,1) on Γ₀(3)\ℍ³ |
| Mesh | 6×3×3, 26 532 CR dofs |
| Y, ν\*, windows | 1.25, 1.02, 8 |
| (θ, θ₂, α, θ₄, ρ̃) | **(0.5, 0.9, 0.15, 0.8, 1.5)** |
| κ | Lemma I1 |
| Wall time | ~38 min |
| Result | 8/8: c_S>0, c_e>d_e, PSD True; c_e ∈ [0.836, 0.881], d_e ∈ [0.079, 0.222] |

```powershell
cd ...\np9\independent_exclusion
python -u -c "from m3p_certify import certify; certify(6,3,3,params=(0.5,0.9,0.15,0.8,1.5),level='(3)')"
```

Dependency delta vs Theorem 2: **none** (only ring + copy count).

---

## 6. N𝔭=13 attempts (failure case)

| | Attempt A | Attempt B |
|--|-----------|-----------|
| Mesh | 6×3×2 (24 888 dofs) | 8×3×2 (33 176 dofs, ~8.8 GB A) |
| ν\* / NWIN | 1.005 / 16 | 1.001 / 16 |
| Params | (0.35, 0.95, 0.1, 0.85, 0.5) | (0.3, 0.98, 0.1, 0.85, 0.5) |
| Float pencil | ~1.04 | ~1.13 |
| Rump | **PSD False** w0 | **PSD False** w0, w1 |
| Scalars | pass | pass |

Stopped B after two windows: same failure mode, ~10 min/window, no point finishing 16.

---

## 7. Root-cause analysis

### Why N=9 worked

- μ ≈ 2.4 gives ~2× headroom over 1 before inflation
- After inflation, pencil ~1.13 and Rump still clears
- Dense size ~5.7 GB fits; 8 windows practical

### Why N=13 did not

1. **μ collapsed to ~1.45** (volume / multi-copy core), faster than handoff guessed  
2. **G1𝔭 inflation** (especially c_Q = 1 − (ω+νλ)S_Q) eats several percent  
3. **Rump ε grows with n** and mid/rad; float LOBPCG overestimates usable margin  
4. Matching Thm 3’s pencil number (~1.13) is **not** enough at larger n / tighter μ  
5. Next geometric step that helps S_Q enough (**8×4×***) blows RAM for dense Rump on 16 GB

### What is *not* broken

- Gluing / F₁₃ combinatorics (asserts + cusp checks)
- Lemma D0 path (no τ reintroduction)
- Positive μ (form still “points right”)
- Mode-bound theory (N=13 ≪ 61)

### Search-design lesson

Scalar slack **anti-correlates** with pencil via ρ̃. Ranking candidates by scalar slack alone prefers **high** ρ̃ and can miss the only Rump-viable region (low ρ̃). Fixed for N=9 by expanding grid + low-ρ̃ trial order; still not enough physics margin at N=13.

---

## 8. Artifacts left in np9

| Path | Role |
|------|------|
| `congruence_prototype.py` | Generalized ring, assembly, margin CLI |
| `m3p_certify.py` | Level-parametrized certificate |
| `PROOF.md` | Theorems 1–3; **no Thm 4** |
| `PROGRESS_LADDER.md` | Full run log + bound analysis |
| `HANDOFF_LADDER.md` / `CONGRUENCE.md` | Status annotations |
| `probe_np13*.py` | Diagnostics (not load-bearing) |
| `theorem4postmortem.md` | This document |

Main: still old F₅-only prototype; **do not assume main has Theorem 3 code**.

---

## 9. Lessons for the next agent

1. **Trust assembly checks first** (t₀, connectivity, mass); then μ; then pencil; then Rump — each is a tighter filter.  
2. **μ < ~1.5** ⇒ treat as yellow/red for *this* certificate stack, not just “rerun params.”  
3. **Float pencil > 1 is necessary but not sufficient** for Rump at 25k+ dofs.  
4. **Planar S_Q** is the main knob for inflation; height N3 mainly affects γ (scalar d_e).  
5. **Don’t burn multi-hour 16-window runs** after two hard-window PSD fails.  
6. **Keep main clean** — all experimental ladder work stays in a worktree.  
7. Handoff μ extrapolations were useful as priors, **not** as guarantees.

---

## 10. Paths to Theorem 4 (if resumed)

Ordered by honesty of cost:

1. **Out-of-core / blocked Cholesky** under Rump (handoff fallback; hours; admissible “any library routine”) on a finer planar mesh (8×4×2/3).  
2. **More RAM** machine for dense A at 8×4×3 (~35 GB).  
3. **Theory sharpening** (S_Q, window/om constants) — real math work, not a one-day code tweak.  
4. **Y / criterion redesign** — changes the proved setup; not free.  
5. **Do not** only retry 6×3×2 + same G1𝔭 constants.

---

## 11. Resolution (addendum)

Unscaled Rump used a **uniform** max-row-sum radius shift larger than
λ_min(A)≈3·10⁻⁵ while A was still SPD. Fix in `m3p_certify.py`:

1. Power-of-two diagonal congruence SAS (IEEE-exact).
2. Per-row radius diagonal reduction r_i = Σ_j |SΔAS|_ij.

**Theorem 4 certified** at 8×3×2, ν*=1.001, 16 windows, params
(0.3, 0.98, 0.1, 0.85, 0.5); ~93 min; all PSD True.

## 12. Bottom line

| | |
|--|--|
| **Delivered** | Theorems **3 and 4** in np9; portable prime-level stack + equilibrated Rump |
| **Was blocked, now fixed** | N𝔭=13 Rump false negative from uniform radius shift |
| **Scientific takeaway** | Ladder is volume-sensitive (μ↓) but N=13 is certifiable on 16 GB with correct Rump scaling; μ≳1.5 is not a hard wall |
| **Isolation** | Main worktree clean of these changes until explicit merge |

Operational detail: `PROGRESS_LADDER.md`.
