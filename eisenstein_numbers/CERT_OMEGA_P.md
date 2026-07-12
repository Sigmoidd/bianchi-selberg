# Research note: `cert_omega_p` — interval cert for Γ₀(𝔭) over ℤ[ω]

**Status:** **CERTIFIED** Γ₀(N=3) 2026-07-12 (`cert_omega_p.py` 8/8 Rump).  
Paper: `PROOF.md` Theorem B. Pairing matrices: `PAIRING_MATRICES.md`.  
**Typo note:** user “cert_onega_p” → **`cert_omega_p`**.

**Target theorem (rung 1).**  
For \(\mathfrak{p}=(1-\omega)\), \(N(\mathfrak{p})=3\), index \(4\):

> \(\Gamma_0(\mathfrak{p})\backslash\mathbb{H}^3\) has no Laplace eigenvalue in \((0,1)\).

Method: multi-cusp Lax–Phillips + **Theorem G1𝔭** (CR lower bound, Lemma D0)
+ Arb enclosures + Rump BIT 46 — the ω-analogue of `m3p_certify.py`
(Gaussian ladder Thms 2–4), **not** a re-run of Thms 1–4.

---

## 1. Architecture (three layers)

```text
┌─────────────────────────────────────────────────────────────┐
│  Criterion (analytic)                                       │
│  Two cusps ∞,0; |T_∞|=√3/6, |T_0|=3·√3/6  [T0_AREA.md]     │
│  β_α=(1−s)/(|T_α|Y²);  κ_c=1/((1+s)Y²) on t_∞+t_0          │
└───────────────────────────┬─────────────────────────────────┘
                            │
┌───────────────────────────▼─────────────────────────────────┐
│  Theorem G1𝔭 (discrete)                                     │
│  Reference-Cell Principle: one P₃ cell → NC=4 isometric     │
│  copies; conforming glue on cross-copy faces; Neumann on     │
│  self-IDs.  Lemma D0: σ_h = √(index)·α_h^ref  (no τ in d_e)  │
│  N_h − D_h ⪰ 0, c_e ≥ d_e on windows covering (0,1)          │
└───────────────────────────┬─────────────────────────────────┘
                            │
┌───────────────────────────▼─────────────────────────────────┐
│  Machine (cert_omega_p)                                     │
│  Arb mid/rad on ref cell → exact scatter → Rump SAS/row     │
└─────────────────────────────────────────────────────────────┘
```

Parallel documents:

| Layer | Gaussian | Eisenstein ω |
|-------|----------|----------------|
| Criterion | `CONGRUENCE.md` §2, `0cuspchecks` | `T0_AREA.md`, `cusps_omega.py` |
| G1𝔭 | `CONGRUENCE.md` §7 | **to state** (this note §3) |
| Float multi-copy | `congruence_prototype.py` | `congruence_omega_proto.py` |
| Interval cert | `m3p_certify.py` | **`cert_omega_p.py` (planned)** |
| Level-1 only | `m3_certify.py` | `cert_omega.py` (**done**) |

---

## 2. Component inventory (reuse vs rewrite)

### 2.1 Reuse as-is (from independent_exclusion)

| Piece | Source | Role |
|-------|--------|------|
| `tet_arb_data`, `mid_rad`, `upper`, `amax/amin`, `a_yf` | `m3_certify` | tet geometry in arb |
| `weighted_moments` | `m3p_certify` | exact-weight Q/M/a enclosures |
| `rump_certify_inplace` + `power_of_two_diag_scale` | `m3p_certify` | SAS equilibration + per-row radii |
| Rump BIT 46 theory | same | PSD certificate |

### 2.2 Reuse with |T| / index adaptation (from eisenstein_numbers)

| Piece | Source | Change for cert_omega_p |
|-------|--------|-------------------------|
| `build_P3_cell` | `reference_cell` | mesh freeze N_tri=6, N3=3 |
| `orient_tets`, Lemma G arb | `cert_omega` | same, on ref cell |
| Cusp areas / β | `cusps_omega`, `T0_AREA` | proved; freeze in code |
| Residue P¹ + gluing perms | `residue_omega` | exact combinatorics |
| Face→face maps | `face_pairings_p3` | must be cert-frozen (see §5) |
| Float assembly pattern | `congruence_omega_proto` | lift mid/rad scatter |

### 2.3 Rewrite / new (G1𝔭 ω form)

| Piece | Why new |
|-------|---------|
| `window_coeffs_p` | β_∞, β_0 use \|T_α\| not Gaussian 2/Y² and 2/(N Y²) with \|T\|=1/2 |
| `sigma_h` | √4 · α_h^ref (index 4), not √6 / √(N+1) for Gaussian N=5 |
| `build_A` two-cusp | D has β̃_∞ t_∞t_∞ᵀ + β̃_0 t_0t_0ᵀ (rank-two boundary) |
| `scaled_radius_rows` | two t-vectors (already in m3p — port) |
| τ in level-1 cert | **omit from d_e** if Lemma D0 adopted (m3p path) |

---

## 3. Theorem G1𝔭 for ω (coefficients to implement)

Match `CONGRUENCE.md` §7 / `m3p_certify.window_coeffs`, with:

```text
λ⁺     = 1 − (s⁻)²
β_∞⁺   = (1 − s⁻) / (|T_∞| Y²)     # |T_∞| = √3/6
β_0⁺   = (1 − s⁻) / (|T_0|  Y²)     # |T_0|  = N |T_∞|
Δκ     = max_W |κ_c(s) − κ_c(s₀)|,  κ_c = 1/((1+s)Y²)
d₂     = ρ̃ / θ₂,   ρ̃ = ρ(1−θ₂)
ω      = 12 d₂ V_S                   # same Young bookkeeping as G1p
c_Q    = 1 − (ω + ν* λ⁺) S_Q
c_Σ    = 1 − (ω + ν* λ⁺) S_S
λ̃      = ν* λ⁺ (1+S_M) + ω S_M
β̃_∞    = ν* β_∞⁺ + 4 d₂ Δκ²
β̃_0    = ν* β_0⁺ + 4 d₂ Δκ²
σ_h    = √(N(𝔭)+1) · α_h^ref         # Lemma D0; index=4 for N=3
c_e    = c_Q · (1 − (1/θ₄−1) ρ_w) − ρ̃(1/θ−1) σ_h²
d_e    = λ̃ (1+1/α) γ²                # no τ term under D0

N_h = c_Q [Q_pc + (1−θ₄) Q_rem] + ρ̃(1−θ) ℓ₀ ℓ₀ᵀ
D_h = λ̃(1+α) M_ex + β̃_∞ t_∞t_∞ᵀ + β̃_0 t_0t_0ᵀ
ℓ₀  = a + κ_c(s₀)(t_∞ + t_0)
```

**Lemma D / D0 — checked for this setup (see §3.1).** Holds for P₃ multi-copy
tops; adopt D0 (drop τ from d_e). Level-1 `cert_omega` still uses classical
G1 with τ (valid, lossier).

---

## 3.1 Lemma D / D0 in the ω · P₃ · multi-copy setting

Source statement (`CONGRUENCE.md` §7, abbreviated):

> **Lemma D.** Cusp cross-sections at height \(Y\) are exactly tiled by
> boundary faces of the mesh, and \(I_{\mathrm{CR}}\) matches the mean of
> \(v\) on every mesh face. With \(e=v-I_{\mathrm{CR}}v\),
>
> - **(D0)** \(t_\infty(e)=t_0(e)=0\), i.e. \(t_\alpha(I_{\mathrm{CR}}v)=t_\alpha(v)\) exactly;
> - **(D3)** \(|a(e)|\le \sqrt{N_C}\,\alpha_h^{\mathrm{ref}}\sqrt{Q_{\mathrm{pc}}(e)}\)
>   (disjoint copies),
>
> so \(\sigma_h=\sqrt{N_C}\,\alpha_h^{\mathrm{ref}}\) and \(d_e\) has **no**
> trace/\(\tau\) part.

Also the level-1 remark (`lower_bound_theory.md`): if \(\{y=Y\}\) is tiled by
mesh faces and \(I_{\mathrm{CR}}\) has zero face-mean error on those faces,
then \(t(e)=0\) and \(\tau_h=0\).

### Hypotheses vs our mesh (checked)

| Hypothesis | Gaussian box | ω P₃ multi-copy | Verdict |
|------------|--------------|-----------------|--------|
| Top plane \(\{y=Y\}\) is a mesh boundary | yes (flat top) | **yes** — all top nodes at \(y=Y\) exactly | OK |
| Top tiled by CR faces (triangles) | yes | **yes** — `top_faces` Kuhn/prism tris; \(\sum\mathrm{area}=\|T\|\) | OK |
| \(I_{\mathrm{CR}}\) matches face means | CR definition | **same CR space** on tets | OK |
| \(t_\alpha=\sum\) top integrals over class-\(\alpha\) copies | Facts A–B | **same** (T0_AREA chimneys; 1+3 copies) | OK |
| Tops **not** identified by side gluing | box pairings exclude top | **yes** — T1/Tw/U/S pair maps never include TOP fids (checked in code) | OK |
| Copies disjoint interiors | yes | yes (glue only boundary faces) | OK |
| Floor / sphere in (D0)? | no — only tops | no — S pairs floor only | N/A |

**Numerical spot-check (mesh 6×3):**  
`top y min=max=Y=1.25`; sum of top triangle areas \(=\sqrt3/6\); pair maps hit 0 TOP faces.

### Conclusion for `cert_omega_p`

| Claim | Status in this context |
|-------|------------------------|
| **(D0)** \(t_\alpha(I_{\mathrm{CR}}v)=t_\alpha(v)\) | **Holds** under the standard CR face-mean definition + horizontal top tiling. Independent of P₃ vs rectangle **as long as tops stay free faces at \(y=Y\)**. |
| **(D3)** / \(\sigma_h=\sqrt{N_C}\alpha_h^{\mathrm{ref}}\) | **Holds** with \(N_C=N(\mathfrak{p})+1=4\) (same disjoint-copy Cauchy–Schwarz as √6 for Gaussian N=5). |
| Drop \(\tau\) from \(d_e\) | **Correct for G1𝔭**; do **not** reintroduce slab τ (HANDOFF_LADDER warning). |
| Level-1 `cert_omega` | Still uses τ>0; valid a fortiori; **not** the G1𝔭 path. |

### What Lemma D does *not* buy

- Does **not** justify vertical face pairings (that is Poincaré / `face_pairings_p3`).
- Does **not** replace Lemma G (floor inclusion) or Lemma S (sliver).
- Does **not** remove \(\alpha_h\) from \(\sigma_h\) — volume functional error \(a(e)\) remains; only the **trace** part of \(\mathcal{L}\) vanishes.
- Gluing must remain conforming on side faces; broken glue can disconnect the mesh but does not by itself falsify (D0) on each copy’s top.

### Implementation rule

In `window_coeffs_p` / Rump pencil for Γ₀:

```text
d_e = λ̃(1+1/α) γ²          # NO + β̃ τ² terms
σ_h = √(index) · α_h^ref    # index=4 for N=3
N_h / D_h boundary: uninflated t_∞t_∞ᵀ, t_0t_0ᵀ  (no Young on t)
```

---

## 4. Pipeline steps for `cert_omega_p.certify(q=3, N_tri=6, N3=3)`

| Step | Action | Float analogue | Risk |
|------|--------|----------------|------|
| 0 | Freeze mesh, Y, ρ, θ, θ₂, θ₄, α, ν* | proto defaults | param search may need retune |
| 1 | Build P₃ mesh, orient tets | `build_P3_cell` | low |
| 2 | Lemma G arb on floor faces | `cert_omega.lemma_G_arb` | low (y_f≥√(2/3)) |
| 3 | Ref arb: Q_pc, Q_rem, M_ex, a mid/rad; γ, α_h, ρ_w, S_*, V_S | `m3p.reference_arb` adapted | τ optional if D0 |
| 4 | Top face areas arb → top_pairs | same | low |
| 5 | Face pair maps freeze (dictionary) | `face_pairings_p3` | **medium** — paper lemma |
| 6 | Gluing perms from ℙ¹(𝔽₃) | `residue_omega.gluing_perms` | low (exact) |
| 7 | Union-find global dofs; scatter mid/rad | `congruence_omega_proto` | low |
| 8 | Assert t_∞(1)=|T_∞|, t_0(1)=|T_0| (arb/float) | proto checks | low |
| 9 | Scalar windows c_e≥d_e in arb | new `window_coeffs_p` | medium (σ_h) |
| 10 | Build Ahat; Rump SAS + per-row radii, 8 windows | `m3p` | low at n~5k |

---

## 5. Load-bearing risks (ordered)

### R1 — Face-pairing dictionary (highest math risk)

`face_pairings_p3` pairs RIGHT↔LEFT (T1), LOW↔UP (Tw), etc. by **edge
geometry + height matching**, not a published Poincaré face list with
explicit matrices on each face.

**For a paper cert you need one of:**

1. **Lemma:** those edge pairings are the face pairings of F_3 for the
   generators used in the presentation (cite EGM/Swan + short proof), or  
2. **Computational certificate:** for every paired faces (F,F'), exhibit
   g ∈ Γ with g(F)=F' (matrix with O entries) and check orientation, or  
3. **Relaxation argument:** free BC on unmatched vertical faces only
   *enlarges* H¹ — if still PSD, OK — but then must **not** glue those
   faces (current code glues dictionary matches only).

Float success (μ>0) does **not** replace (1)–(2).

### R2 — Lemma D0 for P₃ tops

**Resolved for this setup (§3.1):** tops are horizontal at \(y=Y\), tiled by
CR faces, not side-glued; (D0) applies. Keep as a cert assert
(`all top nodes y==Y`, `t_α(1)==|T_α|`).

### R3 — σ_h = 2 α_h^ref

Index 4 ⇒ √4=2. Float G1 at level 1 had c_e/d_e ≳ 2; multi-copy inflates
σ_h² by 4 vs one copy — may squeeze c_e. Mitigations: larger ρ̃ carefully,
finer mesh (memory OK), or accept D0 and retune θ.

### R4 — Q_pc / Q_rem split

m3p uses sandwich Q_pc + rem for exact-main-term refinement. Level-1
`cert_omega` used full `I1` stiffness without rem split and still
passed. Prefer **m3p-style split** for G1𝔭 compatibility with published
G1p statement.

### R5 — Gluing exactness

Union-find on (copy, ref_face) must match the mathematical identification
γ_c F_src ∼ γ_j F_dst. Combinatorics = residue action (exact). Geometric
pair_map must match the face pairing isometry’s action on CR face dofs
(face-mean of mapped triangle = paired face). Height-matched pairing is
a **proxy** until R1 is closed.

---

## 6. Memory and time budget (rung 1, mesh 6×3)

| Quantity | Value |
|----------|-------|
| Copies NC | 4 |
| Ref tets | 648 |
| Global CR dofs n | **~5400** |
| Dense float64 A | **~234 MB** |
| Working set (est.) | **~1–2 GB** |
| vs Thm 3/4 | 26–33k dofs, 6–9 GB — **much smaller** |
| 16 GB laptop | **fits with large margin** |
| Arb ref assembly | minutes (648 tets × Taylor) |
| 8× Rump on 5k | seconds–low minutes total |

**Verdict:** memory is **not** the bottleneck. Math/pairing load-bearing is.

---

## 7. Proposed file layout

```text
eisenstein_numbers/
  CERT_OMEGA_P.md          # this research note
  cert_omega_p.py          # implementation (scaffold → full cert)
  congruence_omega_proto.py # float (keep as regression)
  face_pairings_p3.py      # freeze API used by cert
  residue_omega.py
  cusps_omega.py
  T0_AREA.md
```

`cert_omega_p.py` stages (suggested):

1. `reference_arb_p3(N_tri, N3)` — port of m3p.reference_arb to P3  
2. `assemble_global_p(ref, q=3)` — glue + scatter mid/rad, t_inf/t_0  
3. `window_coeffs_p(...)` — §3 formulas  
4. `certify(q=3, ...)` — Lemma G + scalars + Rump loop  
5. CLI: `python -u cert_omega_p.py` / `python -u cert_omega_p.py status`

---

## 8. Implementation phases

| Phase | Deliverable | Exit criterion |
|-------|-------------|----------------|
| **P0** Research | this doc | done |
| **P1** Scaffold | `cert_omega_p.py` checklist + imports + ref_arb stub | runs `status` |
| **P2** Ref arb | Q_pc/rem, M_ex, a, constants on P3 | matches float 1'M1 order |
| **P3** Global scatter | glued mid/rad; t_∞/t_0 exact | asserts areas; n≈5400 |
| **P4** Scalar windows | arb c_e≥d_e all 8 | print table |
| **P5** Rump | 8/8 PSD | **certificate** |
| **P6** Harden R1 | pairing lemma or computational g-matrices | paper-ready |

---

## 9. Open questions for implementation

1. **Adopt Lemma D0?** Recommended yes (align with m3p / smaller d_e).  
2. **Pairing freeze:** dictionary-only vs require explicit g∈SL(2,O) per face pair.  
3. **Default mesh:** 6×3 (float pass) vs coarser first (faster arb). Prefer 6×3.  
4. **Parameter search:** copy m3p float grid over (θ,θ₂,α,θ₄,ρ̃) or freeze from float pencil.  
5. **N=7 next?** Only after N=3 cert; n~10–12k still fine for RAM.

---

## 10. References (in-repo)

- `independent_exclusion/m3p_certify.py` — gold standard multi-copy cert  
- `independent_exclusion/CONGRUENCE.md` §7 — Theorem G1𝔭 + Lemma D0  
- `independent_exclusion/0cuspchecks.md` — cusp-0 audit template  
- `eisenstein_numbers/cert_omega.py` — level-1 Arb+Rump on P3  
- `eisenstein_numbers/T0_AREA.md` — \|T_0\| proof  
- `eisenstein_numbers/congruence_omega_proto.py` — float μ≈1.12 @ 6×3  
- `eisenstein_numbers/face_pairings_p3.py` — edge dictionary  

---

## 11. Bottom line

| Question | Answer |
|----------|--------|
| Is cert_omega_p well-posed? | **Yes** — same shape as m3p + T0_AREA |
| Memory OK? | **Yes** (~0.25 GB dense; ≪ Thm 3/4) |
| Main gap? | G1𝔭 ω coefficients + pairing load-bearing + wiring scatter/Rump |
| First code step? | Scaffold `cert_omega_p.py` Phase P1–P2 (ref arb on P3) |
