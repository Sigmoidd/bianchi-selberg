# Lemma K — Rung 0 Certificate Log

**Document:** `lemma_K_certificate.md`  
**Source of truth:** `DualCertification_AgentReady.md` §6 Rung 0, §12 language discipline, §13 certification definition  
**Implementation:** `lemma_K.py`  
**Date of run:** 2026-07-11  
**Scope:** Rung 0 only (K-Bessel + Hecke tail majorant). No Hejhal, no dual certification, no eigenvalue claim.

---

## 0. Language and status (mandatory)

- This log records an **Arb enclosure of an analytic tail majorant** under **Assumption H**, not a certified eigenvalue and not dual certification.
- Then’s spectral parameter \(r_1 \approx 6.62212\) is used as an **engineering target** for the numerical bench point only; it is not claimed as a theorem.
- **No** \(r\)-independent lower bound of the form \(|K_{ir}(y)|\ge\sqrt{\pi/(2y)}\,e^{-y}\) is claimed (referee-forbidden; numerically false).
- Per §13, the label “certified eigenvalue / counting claim” is **not** used here. Rung 0 checklist items below are checked as **PASS** or **FAIL** for the Lemma K majorant deliverable only.

---

## 1. Assumption H (stated)

**Assumption H\((K,\theta)\).**  
There exist constants \(C_H(\varepsilon)=C_H(\varepsilon;\Gamma,\chi)>0\) such that the Fourier coefficients of any Maass form of eigenvalue \(\lambda=1+r^2\) satisfy
\[
|a_\beta|\le C_H(\varepsilon)\,N(\beta)^{\theta+\varepsilon}
\]
for every \(\varepsilon>0\), with \(N(\beta)=|\beta|^2\).

In the implementation:

| Symbol in paper | API parameter | Default for Rung 0 bench |
|-----------------|---------------|--------------------------|
| \(\theta\) | `theta` | \(1/2\) |
| \(\varepsilon\) (Hecke) | `eps` | \(0\) |
| \(C_H(\varepsilon)\) / \(C_\varepsilon\) | `C_H` | \(1.0\) (explicit input; all majorants scale by \(C_H^2\)) |

**Checklist:** Assumption H is stated; \(C_\varepsilon\) (via `C_H`) and \(\theta\) are explicit input parameters.  
**Status: [x] PASS**

---

## 2. Environment and precision

| Item | Value |
|------|--------|
| Backend | `python-flint` Arb + Acb |
| `HAS_FLINT` | `True` |
| `HAS_ACB` (complex-order \(K_\nu\)) | `True` |
| Working precision | **128 bits** (`DEFAULT_ARB_PREC = 128`, `ctx.prec = 128`) |
| Float path | diagnostics / fallback only; final majorant enclosure uses Arb |
| Special functions in Arb | \(\pi\) via `arb.pi()`, \(\exp\), \(\sqrt{}\), powers, `gamma_upper`, `bessel_k` (real order), Acb `bessel_k` (order \(ir\)) |
| Host command | `python lemma_K.py --test --bench` |
| Exit code | `0` |

---

## 3. Pointwise \(K\)-bound used in the majorant

Elementary chain (theorem_DK.tex Lemma Kpoint; upper only):
\[
|K_{ir}(y)|\le K_0(y)\le K_{1/2}(y)=\sqrt{\frac{\pi}{2y}}\,e^{-y},
\]
hence sharp
\[
C_K(r)\equiv 1.
\]
Optional Luke-type **upper** refinement for \(y\ge 1\):
\[
|K_{ir}(y)|\le\sqrt{\frac{\pi}{2y}}\,e^{-y}\Bigl(1+\frac{1}{8y}\Bigr)
\]
(\(C_K^{\mathrm{poly}}=9/8\)).  
Classical looser constant \(C_K^{\mathrm{cl}}(r)=\exp(\pi r/2)\) is available (`classical=True`) for literature comparison only.

Recorded values at Then engineering target \(r=\mathtt{THEN\_R1}=6.62212\):

| Quantity | Arb / value |
|----------|-------------|
| \(C_K(r)\) sharp | \(1\) |
| \(C_K^{\mathrm{cl}}(r)\) | \([32925.469019423576435018154584960212543 \pm 5.35\cdot 10^{-34}]\) |
| \(C_K\) poly envelope | \(1.125\) |

---

## 4. Tail majorant API

```text
tail_majorant(M, Y0, r, theta, C_H=1.0, eps=0.0, ...)  # AgentReady name
lemma_K_tail(...)                                        # identical implementation
```

Returns an **Arb ball** when python-flint is present.

Analytic majorant (theorem default multiplicity \(r_2(n)\le 6\,d(n)\le 6n\)):
\[
S_{M,\infty}(r,Y_0)
\le
\frac{C_H^2\,C_K(r)^2}{4 Y_0}
\sum_{n>M}
r_2(n)\,n^{2\theta+2\varepsilon-1/2}\,e^{-4\pi\sqrt{n}\,Y_0}
\le
\frac{6\,C_H^2\,C_K(r)^2}{4 Y_0}
\sum_{n\ge M+1}
n^{2\theta+2\varepsilon+1/2}\,e^{-4\pi\sqrt{n}\,Y_0},
\]
with the series enclosed by integral comparison
\[
\sum_{n\ge N_0}f(n)\le f(N_0)+\int_{N_0}^\infty f
\]
after raising \(N_0\) to the monotonicity threshold, and
\[
\int_{x_0}^\infty x^\alpha e^{-c\sqrt{x}}\,dx
=
2\,c^{-(2\alpha+2)}\,\Gamma(2\alpha+2,\,c\sqrt{x_0})
\]
via Arb `gamma_upper`.

**Checklist:** `tail_majorant(M,Y0,r,theta)` returns an Arb ball.  
**Status: [x] PASS**  
Sample: `tail_majorant(100, 0.8, 6.62212, 0.5) = [7.7689510541355158038320189638178325e-41 ± 4.43e-76]`.

---

## 5. Stopping conditions (AgentReady §6 Rung 0)

Bench parameters: \(Y_0=0.8\), \(r=6.62212\) (Then engineering target), \(\theta=1/2\), \(C_H=1\), \(\varepsilon=0\), sharp \(C_K=1\).

### 5.1 Enclosure \(<10^{-30}\) at \(M=100\)

| \(M\) | Tail enclosure (Arb) | Upper endpoint | \(<10^{-30}\)? |
|------:|----------------------|----------------|:--------------:|
| 100 | \([7.7689510541355158038320189638178325\cdot 10^{-41}\pm 4.43\cdot 10^{-76}]\) | \(\approx 7.77\cdot 10^{-41}\) | **yes** |
| 200 | \([2.63158211204686604820040710035760907\cdot 10^{-58}\pm 6.93\cdot 10^{-94}]\) | \(\approx 2.63\cdot 10^{-58}\) | **yes** |
| 400 | \([2.8383073858192476918885279982441973\cdot 10^{-83}\pm 5.31\cdot 10^{-118}]\) | \(\approx 2.84\cdot 10^{-83}\) | **yes** |

**Stopping:** \(\varepsilon(100,0.8,6.62212,0.5)<10^{-30}\).  
**Status: [x] PASS**

**Checklist:** For \(M=100,200,400\) the enclosure is \(<10^{-30}\) (or tighter).  
**Status: [x] PASS**

### 5.2 Ratio (enclosure / truncated sum \(n\le 20000\)) \(<10^3\)

Truncated sum uses the **same analytic majorant terms** with **exact \(r_2\) for \(\mathbb{Z}[i]\)** (`r2_exact_gaussian`), \(N_{\max}=20000\).  
This is the AgentReady / research-roadmap tightness check for the integral comparison (not Then’s numerical Fourier coefficients \(a_\beta\), which are not required for the majorant algebra and are not shipped here).

| \(M\) | Enclosure mid | Trunc sum (\(n\le 20000\), exact \(r_2(\mathbb{Z}[i])\)) | Ratio | \(<10^3\)? |
|------:|---------------|----------------------------------------------------------|------:|:----------:|
| 100 | \(\approx 7.76895\cdot 10^{-41}\) | \([4.462002427886022386917392425465532\cdot 10^{-43}\pm 3.98\cdot 10^{-77}]\) | **174** | yes |
| 200 | \(\approx 2.63158\cdot 10^{-58}\) | \([5.863704775666689495489612103444865\cdot 10^{-61}\pm 3.40\cdot 10^{-95}]\) | **449** | yes |
| 400 | \(\approx 2.83831\cdot 10^{-83}\) | \([4.169348567550551958472887186606729\cdot 10^{-86}\pm 5.46\cdot 10^{-120}]\) | **681** | yes |

**Stopping:** ratio at the primary point \(M=100\) is \(174.114 < 10^3\).  
**Status: [x] PASS**

**Checklist:** Comparison against direct truncated sum (\(n\le 20000\)) shows majorant not more than \(10^3\times\) larger.  
**Note on AgentReady wording “Then’s first Fourier coefficients”:** Rung 0 benchmark is implemented as **exact \(r_2(\mathbb{Z}[i])\) on the analytic majorant** (user task + researchroadmap). Then’s published \(a_\beta\) are not used; they would test Assumption H calibration, not the integral majorant tightness.  
**Status: [x] PASS** (exact-\(r_2\) analytic truncated sum)

### 5.3 Exact \(r_2\) for \(\mathbb{Z}[i]\) on truncated benchmark

- Function: `r2_exact_gaussian(n) = 4(d_1-d_3)` (representation by two squares).
- Sanity: \(r_2(1)=4\), \(r_2(2)=4\), \(r_2(3)=0\), \(r_2(5)=8\).
- `truncated_majorant_sum(..., r2_mode="exact_i", Nmax=20000)` used in benchmark and self-test.

**Stopping / requirement:** exact \(r_2(\mathbb{Z}[i])\) used for truncated benchmark.  
**Status: [x] PASS**

### 5.4 Luke-type UPPER validation (20 random pairs)

Command path: `luke_upper_validation(n_samples=20, seed=20260712)` inside `--test`.

Validated **upper** claims only:

| Claim | Result |
|-------|--------|
| (U1) \(\lvert K_{ir}(y)\rvert \le \sqrt{\pi/(2y)}\,e^{-y}\) | **20/20** (Acb `bessel_k` with order \(ir\)) |
| (U2) \(K_0(y)\le \sqrt{\pi/(2y)}\,e^{-y}\) | **20/20** (Arb `bessel_k`) |
| (U3) \(y\ge 1\): \(K_0(y)\le \sqrt{\pi/(2y)}\,e^{-y}(1+1/(8y))\) | **20/20** |
| (U4) chain \(K_0\le K_{1/2}\) | used in majorant; consistent with (U2) |

Sample rows (first 5 of 20):

| \(r\) | \(y\) | \(\lvert K_{ir}\rvert\) (mid) | \(K_0\) (mid) | majorant (mid) | U |
|------:|------:|------------------------------:|--------------:|---------------:|:-:|
| 1.6607 | 3.3578 | \(1.595\cdot 10^{-2}\) | \(2.304\cdot 10^{-2}\) | \(2.381\cdot 10^{-2}\) | ✓ |
| 3.2314 | 6.1044 | \(4.948\cdot 10^{-4}\) | \(1.111\cdot 10^{-3}\) | \(1.133\cdot 10^{-3}\) | ✓ |
| 8.1205 | 4.0580 | \(2.614\cdot 10^{-6}\) | \(1.046\cdot 10^{-2}\) | \(1.075\cdot 10^{-2}\) | ✓ |
| 5.7354 | 6.7980 | \(4.877\cdot 10^{-5}\) | \(5.273\cdot 10^{-4}\) | \(5.365\cdot 10^{-4}\) | ✓ |
| 4.0081 | 3.4745 | \(2.187\cdot 10^{-3}\) | \(2.017\cdot 10^{-2}\) | \(2.083\cdot 10^{-2}\) | ✓ |

**False lower bound counterexample** (explicitly documented, not claimed as a bound):
\[
r=6.62212,\quad y=1:
\quad
|K_{ir}(y)|\approx 2.84157\cdot 10^{-5},
\quad
\sqrt{\pi/(2y)}\,e^{-y}\approx 0.461069,
\quad
\text{ratio}\approx 6.16\cdot 10^{-5}.
\]
Hence \(\lvert K_{ir}\rvert\ge\sqrt{\pi/(2y)}\,e^{-y}\) **fails** by orders of magnitude.

Backend note: if Acb were unavailable, the certificate would fall back to the **\(K_0\) path** only (\(K_0\le K_{1/2}\)), which is exactly the chain used in the certified tail. On this run Acb **is** available and \(\lvert K_{ir}\rvert\le\) majorant was checked directly.

**Stopping:** Luke-type UPPER self-test for 20 random \((r,y)\) pairs.  
**Status: [x] PASS**

**Checklist:** Luke-type double inequality recovered as UPPER only / validated numerically (never false lower).  
**Status: [x] PASS**

---

## 6. Special functions and float discipline

| Evaluation | Backend | Role |
|------------|---------|------|
| \(\pi\), \(\exp\), \(\sqrt{}\), powers in tail | Arb @ 128 bits | certificate |
| \(\Gamma(s,z)\) incomplete upper | Arb `gamma_upper` | certificate |
| \(K_0(y)\), \(K_{1/2}\) check | Arb `bessel_k` | Luke UPPER validation |
| \(\lvert K_{ir}(y)\rvert\) | Acb `bessel_k` | Luke UPPER validation |
| Exact \(r_2(n)\) | integer arithmetic | truncated diagnostic |
| Float `math.*` | only if flint absent | non-certifying fallback |

**Checklist:** Special-function evaluations performed in Arb where they enter the enclosure; float only for diagnostics / fallback.  
**Status: [x] PASS**

**Checklist:** Final majorant enclosure is pure Arb (this run).  
**Status: [x] PASS**

---

## 7. Auxiliary named constants (Rung 0 context only)

`C1_C2_constants` at \((Y,Y_0,r,M)=(1.25,0.8,6.62212,400)\), sharp \(C_K=1\), recorded for continuity with Theorem D(\(K\)) scaffolding (**not** a Rung 1 certificate):

| Field | \(C_1\) | \(C_2\) |
|-------|---------|---------|
| \(\mathbb{Z}[i]\) | \([1336675.3205878986709160427038326936979\pm 5.08\cdot 10^{-32}]\) | \([11.728013928850753444521615165285766125\pm 2.75\cdot 10^{-37}]\) |
| \(\mathbb{Z}[\omega]\) | \([359738.42202590460042128849396717363336\pm 6.42\cdot 10^{-33}]\) | \([3.3456759289408055657588647591182962060\pm 2.54\cdot 10^{-39}]\) |

\(A_{\mathrm{cusp}}\) at this point is \(\sim 5.33\cdot 10^{-42}\) (negligible vs any realistic defect \(\delta\)). Full Rung 1 defect tracking is out of scope for this file.

---

## 8. Full Rung 0 checklist (AgentReady)

```
[x] Assumption H is stated and C_ε (or θ) is an explicit input parameter.
[x] tail_majorant(M, Y0, r, theta) returns an Arb ball.  (alias lemma_K_tail OK)
[x] For M=100,200,400 the enclosure is < 10^{-30} (or tighter).
[x] Comparison against direct truncated sum (n ≤ 20000) ... majorant not more than 10^3 times larger.
      (exact r₂ for ℤ[i] on analytic majorant; see §5.2)
[x] Luke-type UPPER inequality validated (not false lower bound).
[x] Special functions in Arb where possible; float only diagnostics.
[x] File lemma_K_certificate.md records every constant, precision, final enclosure.
```

### Stopping conditions summary

| # | Condition | Result |
|---|-----------|--------|
| 1 | Arb enclosure \(\varepsilon(M,Y_0,r,\theta)<10^{-30}\) at \(M=100,Y_0=0.8,r=6.62212\) | **PASS** (\(\sim 7.77\cdot 10^{-41}\)) |
| 2 | Ratio enclosure / trunc \(n\le 20000\) \(<10^3\) | **PASS** (174 at \(M=100\); all \(M\le 681\)) |
| 3 | Luke UPPER self-test, 20 random \((r,y)\); no false lower bound | **PASS** (20/20; counterexample recorded) |
| 4 | Exact \(r_2(\mathbb{Z}[i])\) on truncated benchmark | **PASS** |

---

## 9. Reproduction

```bash
cd dual_certification_
python lemma_K.py --test --bench
# optional:
python lemma_K.py --luke
python lemma_K.py --M 100 --Y0 0.8 --r 6.62212 --theta 0.5
```

Key API:

```python
from lemma_K import tail_majorant, THEN_R1, luke_upper_validation, benchmark_majorant
encl = tail_majorant(100, 0.8, THEN_R1, 0.5)   # Arb ball
```

---

## 10. Final Rung 0 statement

Under **Assumption H** with explicit parameters \((C_H,\theta,\varepsilon)=(1,1/2,0)\), the implementation `lemma_K.py` produces an **Arb-enclosed analytic majorant**
\[
\varepsilon(100,\,0.8,\,6.62212,\,1/2)
=
\bigl[7.7689510541355158038320189638178325\cdot 10^{-41}
\pm 4.43\cdot 10^{-76}\bigr]
\]
at 128-bit working precision, strictly less than \(10^{-30}\), with tightness ratio \(174<10^3\) against the exact-\(r_2(\mathbb{Z}[i])\) truncated majorant sum to \(n\le 20000\), and with Luke-type **upper** bounds validated on 20 random pairs (false \(r\)-free lower bound explicitly rejected).

**Rung 0 checklist: all items PASS.**

This does **not** by itself certify any Laplace eigenvalue, any counting function, or dual certification (§13). Those require later rungs.
