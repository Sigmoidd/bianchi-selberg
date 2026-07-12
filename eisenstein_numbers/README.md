# Eisenstein numbers — PSL(2, ℤ[ω])

Ring `ℤ[ω]`, `ω = e^{2πi/3}`, and the **Eisenstein–Picard** group
`Γ = PSL(2, ℤ[ω])`.

## Certified level-1 spectral bound (FEM)

**Theorem.** λ₁(PSL(2, ℤ[ω])\ℍ³) ≥ 1  
(no Laplace eigenvalues in (0,1); Neumann core; independent of the trace formula).

| | |
|--|--|
| Writeup | [`PROOF.md`](PROOF.md) |
| Certificate | `python -u cert_omega.py 6 3` |
| Roadmap | [`ROADMAP.md`](ROADMAP.md) |
| Geometry | [`GEOMETRY.md`](GEOMETRY.md), `geometry_fund.py` |

This is **not** Theorems 1–4 (congruence ladder on Γ₀(𝔭) ≤ PSL(2, ℤ[i])).
The ω congruence ladder is still open: [`CONGRUENCE_OMEGA.md`](CONGRUENCE_OMEGA.md).

### Reproduce FEM path

```powershell
cd eisenstein_numbers
python -u geometry_fund.py       # Stage 1 geometry freeze
python -u cr_omega.py             # Stages 5–6 float G1
python -u cert_omega.py 6 3       # Stage 7 Arb + Rump (the cert)
python -u run_roadmap.py          # Stages 1–8 overview
```

### Roadmap status

```
Geometry → Reference cell → Float M0 → Constants
  → CR → Window theorem → Interval cert → Congruence ladder
```

| Stage | Status |
|-------|--------|
| 1–2 Geometry + P₃ cell | **done** |
| 3–4 M0 + constants | **done** (float) |
| 5–6 CR + float G1 | **8/8 PASS** |
| 7 Interval cert | **8/8 CERTIFIED** (`cert_omega.py`) |
| 8 Congruence Γ₀(𝔫) | scaffold only |

Frozen cert mesh: P₃, `N_tri=6`, `N3=3`, `Y=1.25`, `ρ=55`, 8 s-windows,
κ = Lemma I1, Rump with SAS equilibration + per-row radii.

---

## Relation to the rest of the repo

| Path | Role |
|------|------|
| `../independent_exclusion/` | Picard level-1 + **Thms 1–4** (ℤ[i] congruence) |
| `../verify_eisenstein.py` | Mechanical constants (vol, systole, L(1,χ₋₃), …) |
| `../bianchi_omega*.py`, `../cuspidal_ce.py` | STF / CE engine for ω |
| **this directory** | FEM independent exclusion for ω + ring/CE sandbox |

| Path | Spectral method | Claim |
|------|-----------------|--------|
| `framework.py` + STF | Trace formula / CE | software gate B&lt;1 (see `AUDIT.md`) |
| `cert_omega.py` | FEM + Rump | **machine theorem** λ₁ ≥ 1 (this cert) |

---

## STF / CE framework (separate from FEM)

```powershell
python -u smoke_test.py
python -u framework.py          # CE extension gates
```

| Step | What |
|------|------|
| 0 | ℤ[ω] ring smoke |
| 1 | Picard CE gate (must pass first) |
| 2–5 | Eisenstein CE enumeration + Arb B |

Last known: B(Q(ω)) ∈ [0.525, 0.544] &lt; 1. Load-bearing math for paper CE:
[`AUDIT.md`](AUDIT.md).

---

## Theorems (FEM path)

| Theorem | Claim | Reproduce |
|---------|--------|-----------|
| **A** | λ₁(PSL(2,ℤ[ω])\ℍ³) ≥ 1 | `python -u cert_omega.py 6 3` |
| **B** | no eigenvalue in (0,1) on Γ₀(1−ω)\ℍ³ | `python -u cert_omega_p.py 3 6 3` |

Writeup: [`PROOF.md`](PROOF.md). Pairing matrices: [`PAIRING_MATRICES.md`](PAIRING_MATRICES.md)
(`python -u pairing_matrices.py`). Live checklist: [`PROGRESS.md`](PROGRESS.md).

## Optionals (done)

| Item | Artifact |
|------|----------|
| Non-Neumann FE + κ=CZZ | `python -u non_neumann_omega.py 6 3` → 8/8 float |
| Journal draft | `papers/paper3_eisenstein.tex` |

## Remaining work

1. Ladder rungs N=7, 13 (float then cert).

---

## Notes

- Prefer Arb (`python-flint`) for anything certificate-bound.
- Default FEM domain is **EGM P₃** (`domain="P3"`). Legacy R_comp
  parallelogram is M0-only (degenerates at |z|=1).
- Float G1: `cr_omega.py`. Interval cert: `cert_omega.py`.
