# Dual certification — chunked execution plan

**Source of truth:** `DualCertification_AgentReady.md`  
**Rule:** finish one chunk’s checklist before starting the next. Never label “certified” unless §13 holds.

---

## Chunk map

```text
CHUNK 0 ── Rung 0: Lemma K
CHUNK 1 ── Rung 1: Theorem D(K) + defect_bound_arb
CHUNK 2 ── Milestone 2: κ vs M (+ preconditioner)
CHUNK 3 ── Route A counting prototype (level 1)
CHUNK 4 ── Single-cusp interval Hejhal + D(K) plug-in
CHUNK 5 ── Rung 2 dual overlap (FEM ∩ Hejhal) + language
CHUNK 6 ── Route B / counting certificate
CHUNK 7 ── Rung 3 two-cusp coupling + Krawczyk N=5
CHUNK 8 ── Rung 4 N=5 dual cert
CHUNK 9 ── Rung 5 width-vs-N + Rung 6 field port
```

Dependencies: `0 → 1 → {2,3,4}`; `1+3+4 → 5 → 6`; `1 → 7 → 8`; `8 → 9`.

---

## CHUNK 0 — Lemma K (AgentReady Rung 0)

| | |
|--|--|
| **Unlocks** | All defect bounds |
| **Deliverables** | `lemma_K.py`, `lemma_K_certificate.md` |
| **Stop when** | All green below |

**Checklist**
- [ ] Assumption H stated; `θ` / `C_H` explicit inputs
- [ ] `tail_majorant(M,Y0,r,theta)` → Arb ball
- [ ] M=100,200,400 at Y0=0.8, r≈6.622: enclosure `< 1e-30`
- [ ] Enclosure / trunc sum (n≤20k, exact r₂ Z[i]) ratio `< 1e3`
- [ ] Luke **upper** only validated (no r-free lower bound)
- [ ] Final numbers in Arb; float diagnostics only
- [ ] `lemma_K_certificate.md` records constants + precision

**Commands**
```bash
cd dual_certification_
python lemma_K.py --test --bench
```

---

## CHUNK 1 — Theorem D(K) (AgentReady Rung 1)

| | |
|--|--|
| **Depends on** | Chunk 0 |
| **Deliverables** | `theorem_DK.tex`, `constants_DK.md`, `defect_bound_arb.py`, `rung1_certificate.md` |
| **Stop when** | All green below |

**Checklist**
- [ ] Theorem lists Assumptions **H, A, S** + (Y, Y0, r, K)
- [ ] C₁, C₂ Arb-computable named formulas only
- [ ] Trace + Poincaré + y_min metric explicit
- [ ] `|λ̃−λ| ≤ C₁(δ+τ) + C₂ e^{-2π Y0}` plug-in API for Hejhal
- [ ] C₁ ∼ 10⁵–10⁶ at (1.25, 0.8, 6.6); C₂ explicit
- [ ] Tail M=400 `< 1e-70` (negligible)
- [ ] Sensitivity table Y0±0.1, r±0.1 in `constants_DK.md`
- [ ] `defect_bound_arb.py` demo with toy (δ,τ)

**Commands**
```bash
python defect_bound_arb.py --demo
python lemma_K.py --constants
```

---

## CHUNK 2 — Conditioning diagnostic (Milestone 2 / Rung 3 prep)

| | |
|--|--|
| **Depends on** | Chunk 1 (conceptually); can run after 0 |
| **Deliverables** | `hejhal_conditioning.py`, `conditioning_report.md`, CSV/PNG |
| **Stop when** | Table + slope b; **not** full Rung 3 |

**Checklist**
- [ ] log κ(V) vs log M for M=100,200,400,800 at Y0=0.8
- [ ] Diagonal (or block-diagonal) preconditioner D
- [ ] Fit b in log κ ≈ a + b log M **with / without** D
- [ ] Target after preconditioner: **b < 4** (AgentReady stop for Rung 3 analysis)
- [ ] Label float κ as **non-certifying diagnostic**
- [ ] Rung 3 items still open: consistency proof, Krawczyk N=5, two-cusp S radii

---

## CHUNK 3 — Route A counting prototype (§7, queue pri 3)

| | |
|--|--|
| **Depends on** | FEM meshes (existing) |
| **Deliverables** | `route_A_counting.py`, `route_A_status.md` |
| **Stop when** | Coarse smoke run; checklist honest |

**Checklist (from AgentReady Route A)**
- [ ] Truncated core K_Y, certified Y (use Y=1.25)
- [ ] Product-like mesh near artificial boundary **or** error enclosed
- [ ] Dirichlet + Neumann discrete problems (float first)
- [ ] Bracketing structure λ^N ≤ λ ≤ λ^D documented
- [ ] Artificial-boundary error enclosure — **open until proved**
- [ ] Final N(λ) integer interval — **not certified** until Arb/Rump

**Language:** “addresses the remaining logical gap” — **not** “counting certified”.

---

## CHUNK 4 — Single-cusp interval Hejhal (queue pri 4)

| | |
|--|--|
| **Depends on** | Chunks 0–1, ideally 2 |
| **Deliverables** | Hejhal driver + measured (δ,τ) + D(K) application |
| **Stop when** | Interval [c,d] under explicit hypotheses |

**Checklist**
- [ ] M=400, Y0=0.8, Arb ≥128-bit for special functions where used
- [ ] Measured δ, τ recorded
- [ ] `defect_bound_arb` applied → |λ̃−λ| enclosure
- [ ] Preconditioner from Chunk 2
- [ ] Assumption H listed on any existence claim
- [ ] No dual “certified first eigenvalue” without Chunk 5–6

---

## CHUNK 5 — Dual overlap (AgentReady Rung 2 partial)

| | |
|--|--|
| **Depends on** | Chunk 4 + existing FEM exclusion |
| **Deliverables** | Overlap log FEM [a,b] ∩ Hejhal [c,d] |

**Checklist**
- [ ] FEM lower: G1 / CR (0,1) free
- [ ] FEM upper near 45: **engineering target** unless certified Rayleigh
- [ ] Overlap checked in **interval arithmetic**
- [ ] diam(hull) < 1e-4 for level 1 **target**
- [ ] Counting still required for “first eigenvalue” (Chunk 6)

---

## CHUNK 6 — Counting certificate (Route A or B)

| | |
|--|--|
| **Depends on** | Chunk 3 + (optional Hejhal for Route B) |
| **Stop when** | One route’s full checklist passes (§7) |

---

## CHUNK 7–9 — Later (do not start early)

| Chunk | Content |
|-------|---------|
| 7 | Two-cusp well-posedness + Krawczyk N=5 — **DONE** (`rung3_certificate.md`) |
| 8 | N=5 dual pipeline (Rung 4) — **infra DONE**; §13 cert blocked on residual+counting |
| 9 | Width-vs-N (Rung 5) + Z[ω]/inert port (Rung 6) |

---

## Parallelism policy

| Allowed in parallel | Must be serial |
|---------------------|----------------|
| Chunk 2 ‖ Chunk 3 (after Chunk 1) | Chunk 0 before 1 |
| Diagnostics (2) while counting scaffold (3) | Chunk 4 after 1 |
| | Chunk 5 after 4 |
| | Chunk 7 after 1 |

**Current sprint:** finish **Chunk 0 → Chunk 1**, then **Chunk 2 ‖ Chunk 3**.

---

## Agent assignment (this sprint)

| Chunk | Agent focus | Certificate file |
|-------|-------------|------------------|
| 0 | Rung 0 checklist + `lemma_K_certificate.md` | `lemma_K_certificate.md` |
| 1 | `defect_bound_arb.py` + sensitivity | `rung1_certificate.md` |
| 2 | κ vs M preconditioner | `conditioning_report.md` |
| 3 | Route A scaffold | `route_A_status.md` |

Master status rollup: update `STATUS.md` after each chunk closes.

---

## Anti-patterns (agent must refuse)

- Emitting “certified λ₁ = …” without §13 (1)–(8)
- Claiming Luke r-independent lower bound
- Hiding Assumption H
- Calling FEM upper bound near 45 “proved” (use **engineering target**)
- Starting N=5 dual before Chunks 2+7 preconditioning/Krawczyk
