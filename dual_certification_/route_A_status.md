# Route A — Dirichlet–Neumann Bracketing Status

**Source of truth:** [`DualCertification_AgentReady.md`](DualCertification_AgentReady.md) §7 Route A, §11 priority 3.  
**Prototype code:** [`route_A_counting.py`](route_A_counting.py)  
**Date scaffolded:** 2026-07-11  
**Language discipline:** this prototype **addresses the remaining logical gap** of certified counting. It does **not** claim dual certification, a certified first eigenvalue, or \(N(\lambda)=0\).

---

## Headline label

> **engineering diagnostic / addresses counting gap, not certified**

If only float eigenvalues are available (current state): all numerical spectra and candidate brackets below are diagnostics. No Arb/Rump guarantees are claimed.

---

## Route A agent checklist (AgentReady §7)

Copied from the source of truth; marks and notes reflect the **prototype scaffold** only.

```
[x] Compact core K is truncated at a certified height Y with product collar structure.
    note: Mesh infrastructure builds Picard Humbert K_Y at engineering height Y=1.25
    (same Y as independent_exclusion FEM). Product collar above y=Y is the classical
    cusp structure. "Certified" optimal Y / Arb enclosure of geometric constants for
    the *counting* argument is still open — height choice is engineering, reusing
    the exclusion pipeline.

[x] Mesh on K is product-like near the artificial boundary (or the error of
    non-product structure is enclosed).
    note: build_mesh is a tensor product (x1,x2) × height parameter with flat top
    y=Y (product-like collar). Floor is curved y_f with Lemma G lift (polyhedral
    inner approximation). Non-product error near the floor is handled for exclusion
    by lift bounds; for Route A counting we only need product structure near the
    *artificial* boundary, which holds. Full Arb enclosure of any residual
    non-product effect: not done.

[ ] Dirichlet and Neumann eigenvalue problems are both solved with guaranteed
    bounds (CR or Morley + post-processing).
    note: RED / partial. Float D/N + relative mid/rad on Q,M and interval
    residual diagnostics in route_A_arb_scaffold.py. Still no true Arb assembly
    of hyperbolic weights, no Rump, no certified CR GLB. KAPPA1 GLB sketch only.

[ ] The bracketing inequality λ_k^{N,K_h} ≤ λ_k ≤ λ_k^{D,K_h} is proved with
    explicit constants that account for the truncation error.
    note: RED. Candidate float intervals [λ_k^N, λ_k^D] printed. Truncation
    constants scaffolded (truncation_constants_scaffold) but not proved.

[ ] The artificial-boundary error is itself enclosed by an Arb ball (or by a
    comparison with a larger domain).
    note: RED — OPEN. Scaffold collar/strip formulae present. No proved Δ_trunc.

[ ] Final counting function enclosure N(λ) is an integer interval; for the target
    λ it must be [0,0] on (1,λ₁).
    note: YELLOW candidate only. route_A_arb_scaffold produces
    N(λ)∈[# {k:λ_k^D≤λ}, # {k:λ_k^N≤λ}] from float brackets (e.g. N(1)=[0,1]
    includes Neumann null). Not Arb/Rump; not dual-certified.
```

### Green vs red summary

| Item | Status | One-line |
|------|--------|----------|
| Compact core \(K_Y\) + product collar | **GREEN (scaffold)** | Mesh at \(Y=1.25\); product top |
| Product-like mesh near artificial boundary | **GREEN (scaffold)** | Tensor product; flat top face |
| Guaranteed D/N eigenvalue bounds | **RED** | Float only |
| Bracketing + truncation constants | **RED** | Not proved |
| Artificial-boundary error enclosure | **RED / OPEN** | No Arb ball |
| Counting enclosure \(N(\lambda)\in\mathbb{Z}\)-interval | **RED** | Not produced |

---

## Scope of the prototype spectrum

**What is discretized.** Weighted hyperbolic forms on the compact truncated core only:
\[
Q(v)=\int_{K_Y}|\nabla v|^2\,y^{-1}\,dx\,dy,\qquad
M(v)=\int_{K_Y} v^2\,y^{-3}\,dx\,dy.
\]
Artificial BC on the top face \(y=Y\):
- **Neumann:** free CR face means on top;
- **Dirichlet:** top-face CR dofs set to 0.

Side and floor faces use the same Neumann relaxation as the exclusion pipeline (no face pairings). That relaxation is valid for *some* lower-bound strategies when certified; here it is only diagnostic.

**What is not discretized.** Full quotient \(\Gamma\backslash\mathbb{H}^3\), cusp collars \(y>Y\), residual spectrum machinery, Hejhal trial functions, congruence reference-cell gluing.

---

## Relation to Assumption H

Assumption H (Hecke-type coefficient growth) is **not** required for pure variational counting of the Laplacian on a **compact** domain \(K\) with declared BC.

It **is** still relevant when connecting truncated-core counts to the spectrum of the non-compact quotient (Fourier tails / defect bounds). That connection is exactly the open artificial-boundary / truncation item above. Lemma K / Theorem D(K) live in a different part of the dual-certification stack and are not invoked by this scaffold.

---

## Reuse map

| Asset | Role |
|-------|------|
| `../independent_exclusion/cr_prototype.py` | Picard mesh, geometry, CR assembly |
| `../independent_exclusion/m3_certify.py` | Future Arb/Rump patterns (not called yet) |
| `../eisenstein_numbers/reference_cell.py`, `cr_omega.py` | Optional secondary P3 path (`--eisenstein`) |

Primary target: **Picard level-1 \(\mathbb{Z}[i]\)** (AgentReady Rung 2 dependency).

---

## How to run

From this directory:

```bash
# Smoke test (coarse 6×3×3 mesh, <2 minutes)
python route_A_counting.py --smoke

# Slightly denser diagnostic
python route_A_counting.py --N1 8 --N2 4 --N3 4 --neigs 8

# Optional Eisenstein secondary
python route_A_counting.py --smoke --eisenstein
```

---

## Explicit non-claims

1. No dual certification.  
2. No certified \(\lambda_1\) interval near 45 or elsewhere.  
3. No certified \(N(\lambda)=0\) on \((1,\lambda_1)\).  
4. No claim that float CR eigenvalues are rigorous lower/upper bounds.  
5. No Hejhal work (out of scope for this priority-3 scaffold).

When Arb/Rump upgrades land and the artificial-boundary theorem is proved, re-evaluate the checklist items and only then consider a counting certificate under AgentReady §13.
