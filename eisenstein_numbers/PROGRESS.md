# Progress — Eisenstein–Picard FEM

Last updated: **2026-07-12**.

---

## A. Level 1 — PSL(2, ℤ[ω])  ·  **COMPLETE (machine)**

| Item | Evidence |
|------|----------|
| Geometry + P₃ mesh + exact vol | `GEOMETRY.md`, `geometry_fund.py`, `reference_cell.py` |
| Float M0 / CR / G1 | `fem_omega_m0.py`, `cr_omega.py` — **8/8** |
| **Interval cert Arb+Rump** | `cert_omega.py` — **8/8 CERTIFIED** |
| Theorem writeup | `PROOF.md` |

**Claim:** λ₁(PSL(2, ℤ[ω])\ℍ³) ≥ 1 (Neumann core).

```text
python -u cert_omega.py 6 3
```

Paper: Theorems A–B in `PROOF.md` (G4 citations §7). Optional non-Neumann.

---

## B. Stage 8 — congruence ladder Γ₀(𝔫)  ·  **FLOAT RUNG 1 STARTED**

Goal (future): no eigenvalue in (0,1) for Γ₀(𝔭)\ℍ³, 𝔭 ⊂ ℤ[ω].

### Done

| Item | Artifact | Result |
|------|----------|--------|
| Residue fields N=3,7,13 + P¹ gluing | `residue_omega.py` | S²=id, U³=id |
| Combinatorics smoke | `congruence_omega_smoke.py` | PASS |
| Cusp table (∞:0)=(1:N), \|T_∞\|, \|T_0\|=N·\|T_∞\| | `cusps_omega.py` | modes OK @ Y=1.25 |
| **Proof of \|T_0\| formula** | **`T0_AREA.md`** | **proved** (field primes) |
| Multi-copy CR + two-cusp pencil | `congruence_omega_proto.py` | **float** |
| **EGM P_3 face-pairing dictionary** | **`face_pairings_p3.py`** | RIGHT↔LEFT, LOW↔UP, U, S |
| Design doc | `CONGRUENCE_OMEGA.md` | updated |
| **cert_omega_p research** | **`CERT_OMEGA_P.md`**, scaffold `cert_omega_p.py` | architecture + checklist |
| Generator inverse cycles (user ground truth) | `check_generators.py`, `GENERATOR_CYCLES.md` | **PASS**: Γ_∞-eq for \(T_1^{-1},T_\omega^{-1}\); \(S,U\) exact; residue + face dict OK |

### Float evidence — Γ₀(1−ω), N=3, index 4

```text
python -u congruence_omega_proto.py 3 6 3
```

| check | value |
|-------|--------|
| t_∞(1) | = \|T_∞\| = √3/6 ✓ |
| t_0(1) | = \|T_0\| = √3/2 ✓ |
| M-graph components | 1 ✓ |
| min μ (λ∈{0.05,0.5,0.9,0.99}) | **≈ 0.90 > 0** ✓ |

Coarse mesh (4×2): μ fails near λ=1 and/or disconnected.  
Finer mesh (6×3): **PASS float**.

### Interval cert @ 6×3 (`cert_omega_p.py`) — **PASS 2026-07-12**

```text
python -u cert_omega_p.py 3 6 3
# frozen: θ=0.5, θ₂=0.9, α=0.2, θ₄=0.5, ρ̃=9, ν*=1.05
# 8/8 c_e>d_e, 8/8 Rump PSD; n≈5404 dofs, ~234 MB dense
```

| check | value |
|-------|--------|
| Lemma G arb | OK |
| t_∞(1), t_0(1) | exact |
| min c_e/d_e | ≳ 11 |
| Rump PSD | **8/8** |

Pairing matrices: **`PAIRING_MATRICES.md`** + `pairing_matrices.py` **PASS**.

### Proposed rungs

| Rung | 𝔭 | N | index | Status |
|------|---|---|-------|--------|
| 0 | full Γ | 1 | 1 | **certified** |
| 1 | (1−ω) | 3 | 4 | **CERTIFIED** (`cert_omega_p.py` 8/8 Rump) |
| 2 | π\|7 | 7 | 8 | **CERTIFIED** (`cert_omega_p.py 7 6 3` 8/8) |
| 3 | π\|13 | 13 | 14 | combinatorics only |

### Open

| # | Task |
|---|------|
| 1 | ~~Prove cusp areas \|T_α\|~~ **done** — `T0_AREA.md` |
| 2 | ~~EGM wall-pairing dictionary~~ **done** — `face_pairings_p3.py` |
| 3 | ~~cert_omega_p research + implement~~ **CERTIFIED N=3** |
| 4 | ~~Pairing matrices g∈Γ~~ **done** — `PAIRING_MATRICES.md` |
| 5 | ~~Float + cert N=7~~ **done** (float 4×2; cert 6×3) |
| 5b | Float + cert N=13 |
| 6 | 𝔽₄ inert / CRT composites |
| 7 | ~~Non-Neumann FE space + κ=CZZ~~ **done** — `non_neumann_omega.py` 8/8 float |
| 8 | ~~Journal draft~~ **done** — `papers/paper3_eisenstein.tex` |

---

## C. What this is *not*

- Not Theorems 1–4 (ℤ[i] congruence — already done).
- Does not use the Selberg trace formula.
- No interval claim yet for N∈{7,13} (combinatorics only).

---

## D. Next actions (ordered)

1. Ladder: float/cert N=13.
2. Optional: interval cert on the non-Neumann space (implied by Neumann cert).

---

## E. Commands

```text
# Level 1 cert
python -u cert_omega.py 6 3

# Stage 8
python -u residue_omega.py
python -u cusps_omega.py
python -u congruence_omega_smoke.py
python -u face_pairings_p3.py
python -u pairing_matrices.py
python -u check_generators.py
python -u non_neumann_omega.py 6 3           # optional paired FE + κ=CZZ
python -u congruence_omega_proto.py 3 6 3    # float Γ₀(1-w)
python -u cert_omega_p.py 3 6 3              # interval Γ₀(N=3)
```

## F. File map

| File | Role |
|------|------|
| `PROOF.md` | Theorems A (level 1) + B (Γ₀ N=3) |
| `../papers/paper3_eisenstein.tex` | Journal draft (Paper III) |
| `cert_omega.py` / `cert_omega_p.py` | Interval certs |
| `non_neumann_omega.py` | Paired FE space float + κ=I1/CZZ |
| `PAIRING_MATRICES.md` / `pairing_matrices.py` | g∈Γ freeze + checks |
| `cusps_omega.py` / `T0_AREA.md` | Cusp areas + proof |
| `residue_omega.py` | ℙ¹ + gluing perms |
| `face_pairings_p3.py` | EGM edge pairing dictionary |
| `congruence_omega_proto.py` | Multi-copy float CR |
| `CONGRUENCE_OMEGA.md` | Stage 8 design |
