# Rigor closures for the Eisenstein--Picard computation

This note records what is proved, what is cited, and what was changed on
2026-07-09.  It supersedes only the four gaps identified in `HANDOFF.md`.
The trace formula, the standard facts about its scattering determinant, and
the cited finite-order classification remain external mathematical inputs.

## 1. Non-cuspidal elliptic classes

**Status: closed by a published classification, not by the bounded search.**

For \(\Gamma=\operatorname{PSL}_2(\mathbb Z[\omega])\), the finite-order
classification in Elstrodt, Grunewald, and Mennicke, *Groups Acting on
Hyperbolic Space: Harmonic Analysis and Number Theory* (Springer, 1998),
Chapter 4, §4.3, “The singular axes”, gives precisely two singular-axis types
in \(\Gamma\backslash\mathbb H^3\): one order-2 type and one order-3 type.
The order-3 type is cuspidal.  Thus there is exactly one *non-cuspidal*
order-2 conjugacy class.  This is a global statement, so it rules out a class
whose matrices have no small representative.

The role of `elliptic_inventory.py` is therefore verification of the local
data for that class, not a proof of its completeness: it finds the same
order-2 representative, \(|E(R)|=2\), and \(N(T_0)=7+4\sqrt3\).  The NCE
coefficient is consequently

\[
 C_{\rm nce}=\frac{\log(7+4\sqrt3)}8.
\]

This Chapter 4, §4.3 citation must be retained in any standalone proof.  The
in-repository EGM OCR is a related 1982 report and is not used as a substitute
for the global classification theorem.

### Quantitative check

After the tail correction below, the certified upper endpoint is
\(B_\omega<0.543674\), while the one known NCE class contributes
\(0.5094464477\ldots\).  Hence the remaining margin is
\(0.456326\ldots\), or only \(0.8957\ldots\) of one equal NCE contribution.
In particular, the present test function would absorb **zero** additional
classes of the same size: adding one gives the rigorous upper estimate
\(0.543674+0.509447>1.053\).  Positivity alone is therefore not a fallback
proof of NCE completeness.  It is important not to claim otherwise.

## 2. Cuspidal elliptic classes

**Status: closed.**  The following argument turns the finite enumeration into
an exhaustive count.

The field has class number one, so any cuspidal elliptic element is conjugate
to

\[
 g_{u,b}=\begin{pmatrix}u&b\\0&u^{-1}\end{pmatrix},
 \qquad u\in\mathbb Z[\omega]^\times,quad u^2\ne1.
\]

Modulo \(\{\pm I\}\), the two possible nontrivial rotation multipliers are
\(u=\omega,\omega^2\).  Conjugation by a translation replaces \(b\) by
\(b+a(u^{-1}-u)\).  For each of these two \(u\)'s the principal ideal
\((u-u^{-1})\) has norm 3, so every global CE class has a representative in
one of at most \(2\cdot3=6\) residue candidates.

If an element centralizes \(g_{u,b}\), it cannot interchange the two fixed
points: that would conjugate this order-3 rotation to its inverse.  It therefore
fixes infinity and is upper triangular.  Its diagonal entry is a unit \(p\),
and direct matrix multiplication gives the exact criterion

\[
 b\frac{p-p^{-1}}{u-u^{-1}}\in\mathbb Z[\omega].
\]

Checking the six units is a finite exact calculation.  For all six residue
candidates it gives \(|C(g_{u,b})|=6\).  Friedman's Lemma 4.4.4, with
\(k_\infty=l_\infty=1\), \([\Gamma_\infty:\Gamma'_\infty]=3\), and
\(|1-\epsilon^2|^2=3\), now says

\[
 2\sum_i\frac1{3|C(g_i)|}+\frac13=1.
\]

Every actual class has weight \(1/18\), so the identity forces exactly six
actual classes.  Since there were at most six candidates, they are distinct
and exhaustive.  Exact reduction also gives \(|c_i|^2=1\) twice and 3 four
times, hence the CE coefficients \(\frac29\log3\) and \(\frac13\).

## 3. The scattering contour shift

**Status: closed.**  Let

\[
 \Xi_K(s)=s(s-1)|D_K|^{s/2}(2\pi)^{-s}\Gamma(s)\zeta_K(s),\qquad
 C_K=\tfrac12\log|D_K|-\log(2\pi),
\]

and write \(Q=\Xi_K'/\Xi_K\).  The completed zeta is entire, satisfies
\(\Xi_K(s)=\Xi_K(1-s)\), and has no zeros in \(1\leq\Re s\leq2\).
For
\(L(s)=\phi_K'(s)/\phi_K(s)\), substitution of the completion and its
functional equation gives, on \(s=1+ir\),

\[
 L(1+ir)=-Q(1-ir)-Q(1+ir)+\frac1{1-ir}+\frac1{1+ir}.
\]

The \(1/(1-ir)\) term is the pole of
\(\zeta_K(s-1)\) at \(s=2\); this is exactly where that pole is accounted
for.  Since \(h\) is even,

\[
 \frac1{4\pi}\int_{\mathbb R}h(r)L(1+ir)\,dr
 =\frac1{2\pi}\int_{\mathbb R}\frac{h(r)}{1+r^2}\,dr
  -\frac1{2\pi}\int_{\mathbb R}h(r)Q(1+ir)\,dr. \tag{1}
\]

Shift the last integral to \(s=2+it\).  The horizontal sides vanish because
the admissible sinc fourth power is \(O(|t|^{-4})\) in the needed strip and
\(Q=O(\log(2+|t|))\) there.  No zero of \(\Xi_K\) is crossed.  On the new
line,

\[
 Q(2+it)=\frac{\zeta_K'}{\zeta_K}(2+it)
   +\underbrace{\left(\frac1{2+it}+\frac1{1+it}+C_K+\psi(2+it)\right)}_{B_K(t)}.
\]

The absolutely convergent Euler product and the Fourier definition of \(g\)
give

\[
 -\frac1{2\pi}\int h(t-i)\frac{\zeta_K'}{\zeta_K}(2+it)\,dt
   =\sum_{\mathfrak a}\frac{\Lambda_K(\mathfrak a)}{N\mathfrak a}
       g(\log N\mathfrak a).
\]

Thus (1) is precisely the elementary term plus this prime-power term minus
\((2\pi)^{-1}\int h(t-i)B_K(t)dt\), the formula used in
`bianchi_omega_arb.py`.  For \(\mathbb Z[\omega]\), the first ideal norm is
3 and \(\log3>\operatorname{supp}g\), so the prime-power sum is zero.  For
\(\mathbb Z[i]\), the sole surviving term is \((\log2/2)g(\log2)\).

## 4. Explicit tails

**Status: closed after correcting three missing factors of two.**  In
`bianchi_omega_arb.py`, `psi_main`, `pe_main`, and `ps_main` are integrals
over \([-R,R]\).  Their corresponding tail majorants must therefore include
both half-lines.  The previous code included one half-line only.  The three
tail expressions now have the required leading factor 2; the upper endpoint
remains safely below 1.

For \(z=1+ir\), use the first truncation of the digamma asymptotic
expansion (A\&S 6.3.18; DLMF 5.11.2, with the complex remainder bound in
DLMF 5.11(ii)):

\[
 \psi(z)=\log z-\frac1{2z}+\epsilon(z),\qquad
 |\epsilon(z)|\leq\frac{\sec(\tfrac12\arg z)}{12|z|^2}.
\]

Here \(|\arg(1+ir)|<\pi/2\), so the secant is at most \(\sqrt2\).  Also
\(|\log z|\leq\log|z|+|\arg z|\leq\log(1+|r|)+\pi/2\), and \(|z|\geq1\).
Consequently

\[
 |\psi(1+ir)|\leq\log(1+|r|)+\frac\pi2
 +\frac1{2|1+ir|}+\frac{\sqrt2}{12|1+ir|^2}
 <\log(1+|r|)+1+\frac\pi2.
\]

This uses the sector-correct secant factor; the factor-free complex remainder
bound is not valid.  Integration by parts, using
\((\log(1+r))'\leq1/r\), gives the exact closed-form `psi_tail` bracket.

For the shifted line,
\(|\sin(\delta(t-i))|\leq\cosh\delta\), hence
\[
 |h(t-i)|\leq\left(\frac{\cosh\delta}{\delta t}\right)^{2k}.
\]
The recurrence for \(\psi\), the two reciprocal terms, and
\(|C_K|\leq\log(2\pi)\) for \(|D_K|=3,4\) give for \(t\geq R\)

\[
 |B_K(t)|\leq\log(1+t)+1+\frac\pi2+\frac3t+\log(2\pi)
 \leq\log(2+t)+3+\log(2\pi),
\]

provided \(R\ge7\).  The code now rejects smaller \(R\).  The elementary
tail uses \(h(r)/(1+r^2)\ge0\) and
\(\int_R^\infty r^{-2k-2}dr=R^{-2k-1}/(2k+1)\).  These prove every
constant used by the corrected tail bounds.

## Executed checks

On 2026-07-09, `cuspidal_ce.py` was run with the specified local Python 3.13:
it first reproduced the Picard inventory, then gave the Eisenstein multiset
\((|c|,|C|)=\{(1,6)^2,(\sqrt3,6)^4\}\), and finally ran the exact
residue-class/centralizer self-check whose Friedman identity equals 1.
The corrected `bianchi_omega_arb.py` was then run for both fields:

| field | certified enclosure |
|---|---|
| \(\mathbb Q(i)\) | \([0.304576,\ 0.317810]\) |
| \(\mathbb Q(\omega)\) | \([0.525216,\ 0.543674]\) |

Both upper endpoints are below 1.  The wider intervals replace the earlier
one-half-line tail enclosures.
