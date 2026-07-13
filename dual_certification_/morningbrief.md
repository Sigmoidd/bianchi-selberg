# Morning brief — Rung 4 residual path

**Updated:** \(V_h^{P1,\mathrm{per}}\) / \(J_h\) bridge + true-height σ₀ hybrid modes  
**Hard map (do not change):**  
`F5_gluing=true, FEM_exclusion_01=true, two_cusp_hejhal_scan=true, defect_bound_applied=true,`  
`eta_le_eta0=false, width_lt_tol=false, counting_certified=false`  
⇒ **`pipeline_ok=true`, `rung4_certified=false`** ← still correct.

---

## 0. Snapshot

| Shipped | Note |
|---------|------|
| C1 1.0–1.2+**1.6 proved** | C1≈**1.98e3**, η₀≈**6.4e-8**, gap **~3.9** orders @ η=5e-4 |
| Hybrid multi + periodize | best δ≈**0.088** (M=48); larger M worsens absolute δ |
| **`p1_per_jh_bridge.py`** | \(V_h^{P1,\mathrm{per}}\) + \(J_h\) cross-copy defect (**PASS**) |
| **sigma0 hybrid modes** | true-height σ₀ pin; ~10× vs col, not better than periodize |

---

## 1. C1 (closed 1.0–1.2, 1.6)

\(A_{\mathrm{bdry}}=C_{\mathrm{tr}}C_{\mathrm{Sob}}A_{\mathrm{met}}\) (λ-free).  
Open: 1.3–1.5 mesh / directional / Sobolev (~1.2–4× each).

---

## 2. δ_aut residual path

| mode (M=48, jw=2) | δ_aut | fac vs col |
|-------------------|------:|-----------:|
| collocation | 3.39 | 1× |
| multi hybrid | 0.443 | 7.6× |
| **periodize** | **0.088** | **38×** |
| sigma0 (true height) | 0.342 | 9.9× |
| sigma0_per | 0.670 | 5.1× |

True-height σ₀ hybrid is honest production shape but **does not beat** height-matched periodize on absolute δ. Larger M still worsens δ (see `hybrid_scan_M.md`).

---

## 2b. \(V_h^{P1,\mathrm{per}}\) / \(J_h\) bridge (**new**)

**Code:** `p1_per_jh_bridge.py`  
**Level:** Γ₀(2+i), 6 copies, PAIRINGS T1/R/TiR/S from congruence.

| Quantity | Result |
|----------|--------|
| raw CR face dofs (4×2×2 mesh) | 5088 |
| periodic dofs | 4816 (~1.06× reduction) |
| cross-copy merges | 464 edges |
| self-relaxed (Neumann) | 208 |
| \(J_h^{\mathrm{cross}}\)(constant) | 0 |
| \(J_h^{\mathrm{cross}}\)(free Neumann random) | O(1) (~5.9) |
| \(J_h^{\mathrm{cross}}\)(\(V_h^{P1,\mathrm{per}}\) lift) | **0** (by construction) |
| free→periodic L² defect | ~0.24 |

**Semantics (do not violate):**
- Free Neumann Rayleigh is **lower-bound flavor**, not dual upper bound.
- \(V_h^{P1,\mathrm{per}}\) = CR face-means modulo **cross-copy** gluing only.
- Self-faces stay free (Neumann); \(J_h^{\mathrm{self}}\) need not vanish.
- Companion \(J_h^{\mathrm{cross}}(u)\) measures defect of a free trial vs periodic space.
- **Not yet:** Q/M assembly + Rayleigh on \(V_h^{P1,\mathrm{per}}\) for dual cert.

Reproduce: `python p1_per_jh_bridge.py`

---

## 3–4. τ / η₀ / Route A

Unchanged: residual M-scan; U_norm API; Route A trunc RED.

---

## 5. Checklist

### Done
- [x] Hard map honest
- [x] C1 1.1+1.2+1.6 → C1 1.98e3
- [x] Hybrid multi / periodize / hybrid_scan_M
- [x] **\(V_h^{P1,\mathrm{per}}\) identification + \(J_h\) bridge PASS**
- [x] True-height σ₀ hybrid modes (engineering)

### Open
- [ ] Assemble Q/M on \(V_h^{P1,\mathrm{per}}\) + engineering Rayleigh (not dual upper yet without cert)
- [ ] True two-cusp production Hejhal iterate (δ≪0.09)
- [ ] counting_certified / eta_le_eta0 / width_lt_tol / rung4_certified

---

## 6. Next

| Pri | Action |
|----:|--------|
| 1 | Q/M on \(V_h^{P1,\mathrm{per}}\) (reuse congruence assemble_level_p reduction) |
| 2 | Production two-cusp Hejhal residual at eigen-r |
| 3 | Route A Arb D/N |
| 4 | Keep hard map |

---

## 7. Reproduce

```bash
cd dual_certification_
python p1_per_jh_bridge.py
python delta_aut_pairing.py --M 48 --compare --jump-weight 2
python c1_audit_16.py
python validate_gap_numbers.py
python lemma_K.py --test
```

---

## 8. Bottom line

Production C1 **~2e3**, gap **~3.9 orders**. Best Fourier residual still **δ~0.09** (periodize M=48). New **periodic CR identification** is correct and tested: cross-copy \(J_h=0\) on \(V_h^{P1,\mathrm{per}}\), free Neumann has \(J_h^{\mathrm{cross}}=O(1)\). Next engineering step is variational spectrum on that space — still not dual GREEN without residual + counting. Hard map unchanged.

*Path: `dual_certification_/morningbrief.md`*
