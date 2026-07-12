# Generator actions on ideal vertices (checked)

## User data (ground truth)

Action of **inverses** on the Ford/EGM ideal vertices
\(\{\infty,0,1,\omega,\omega^2\}\):

| generator | cycles |
|-----------|--------|
| \(T_1^{-1}\) | \((\infty)(0\ 1)(\omega\ \omega^2)\) |
| \(T_\omega^{-1}\) | \((\infty)(0\ \omega)(1\ \omega^2)\) |
| \(U^{-1}\) (as labeled) | \((\infty)(0)(1\ \omega\ \omega^2)\) |
| \(S^{-1}\) | \((\infty\ 0)(1)(\omega\ \omega^2)\) |

Audit: `python -u check_generators.py` → **OVERALL PASS**.

Always re-check after geometry/generator edits against this table.

## Agreement with our matrices (`geometry_fund.GEN`)

| Generator | Möbius | Match to user cycles |
|-----------|--------|----------------------|
| \(S\) | \(z\mapsto -1/z\) | **Exact** (unit reduction \(-1\sim 1\)): \(S=S^{-1}\) gives \((\infty\ 0)(1)(\omega\ \omega^2)\). |
| \(U=\mathrm{diag}(\omega^2,\omega)\) | \(z\mapsto \omega z\) | User’s labeled \(U^{-1}\) cycle equals our **forward** \(U\): \((1\ \omega\ \omega^2)\). Our matrix \(U^{-1}\) is the inverse cycle \((1\ \omega^2\ \omega)\). Same cyclic group \(\langle U\rangle\); residue gluing already acts by \(\delta^{-1}\). |
| \(T_1,T_\omega\) | \(z\mapsto z+1\), \(z\mapsto z+\omega\) | Pure Möbius leaves the 5-set; images are **Γ_∞-equivalent** to the user section (existence of \(u\in O^*\), \(b\in O\) with \(u\cdot g^{-1}(v)+b=\) user image). Greedy reduction is multi-valued — **user cycles fix the polyhedron section**. Each user section is an involution on the 5-set. |

### What “Γ_∞-equivalent” means

\[
z \sim_{\Gamma_\infty} z' \iff \exists\, u\in O^*,\; b\in O:\quad z = u\, z' + b.
\]

`check_generators.gamma_inf_equivalent` verifies this for every vertex under \(T_1^{-1},T_\omega^{-1},S^{-1},U^{\pm1}\).

## Residue gluing (what the cert uses)

For \(\Gamma_0(\mathfrak{p})\), copies are labeled by \(\mathbb{P}^1(O/\mathfrak{p})\), **not** the 5 vertices.

| Level | Check |
|-------|--------|
| N=3,7,13 | perms bijective; \(S^2=\mathrm{id}\); \(U^3=\mathrm{id}\); round-trip \(p\cdot\delta^{-1}\cdot\delta=p\) |
| N=3 special | \(\omega\equiv 1\) in \(\mathbb{F}_3\), so \(T_1\equiv T_\omega\) as residue translations; \(U\equiv\mathrm{id}\) on \(\mathbb{P}^1(\mathbb{F}_3)\). Gluing cross-edges come from \(T_1/T_\omega\) and \(S\) only — **expected**. |

## Face dictionary vs vertex cycles

`face_pairings_p3` pairs **edges of \(P_3\)** (RIGHT↔LEFT, LOW↔UP, …).  
Vertex cycles constrain how generators move **ideal vertices** of \(F_3\).  
Both are compatible with the same generators; the cert’s machine glue is
residue \(\mathbb{P}^1\) + edge dictionary (Lemma D0 tops free).

## Cert status

Machine cert for \(\Gamma_0(1-\omega)\) stands (`cert_omega_p.py` 8/8).  
Pairing matrices: `PAIRING_MATRICES.md` / `pairing_matrices.py` (**PASS**).  
Paper theorems: `PROOF.md` Theorems A–B.
