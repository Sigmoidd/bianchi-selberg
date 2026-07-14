# Track B projected-mass certificate

**Status:** the Track B projected-mass input is interval certified using an
exact cusp plateau. The full Theorem D(K) residual is not yet certified, so
`rung4_certified=false` remains mandatory.

## Certified conclusion

For the fixed `M=8`, `r=6.62208` six-copy finite trial, define

\[
W_B=P_{Q,-}P_{-z,+}P_{\rm old}F,
\qquad
P_{Q,-}=\frac{I-Q}{2},
\qquad
P_{-z,+}f=\frac{f(z)+f(-z)}2.
\]

On the Picard cusp slab

\[
\Omega_B=[-1/2,1/2]\times[0,1/2]\times[1.20,1.45],
\]

the canonical Arb run proves

\[
\boxed{\|W_B\|_{L^2(\Omega_B)}>0.20380681000849094.}
\]

The exact recorded coefficient decimals are rigorously enclosed by Arb. The
source trial SHA-256 is
`c77b711b500cae0442376eec1f9ba79f0fb66889e6c3fdd67dac6c6c5e7c733a`.

## Exact six-copy collapse

For a zero-cusp mode `(a,b)`, the five translated components have frequency

\[
(u,v)=((2a-b)/5,(a+2b)/5).
\]

Their translation sum vanishes unless both `u` and `v` are integers. The
surviving modes are in bijection with the 24 infinity-cusp modes. For an
infinity mode `(u,v)`, the unique lift is

\[
(a,b)=(2u+v,-u+2v),
\qquad a^2+b^2=5(u^2+v^2).
\]

Hence the oldspace coefficient is computed without sampled cancellation:

\[
d_{u,v}=\frac{c^{\infty}_{u,v}+5c^0_{2u+v,-u+2v}}6.
\]

The script then forms

\[
e_{u,v}=\frac12\left(
 \frac{d_{u,v}+d_{-u,-v}}2
 -\frac{d_{-v,u}+d_{v,-u}}2
\right).
\]

This gives exact `-z` evenness and quarter-turn oddness up to outward Arb
rounding balls of radius below `4.0e-53` at 192 bits.

## Two independent norm computations

| Bits | Height segments | Parseval/evenness lower |
|---:|---:|---:|
| 128 | 64 | 0.0164629141 |
| 160 | 128 | 0.1295707099 |
| 192 | 256 | 0.1794699173 |
| 192 | 512 | **0.2038068100** |

The first method uses exact full-torus Parseval and divides by two because
`W_B(-z)=W_B(z)`. The second directly integrates the half-period Fourier Gram
and certifies every block eigenvalue using Rump's method. Its lower bound is
weaker because it replaces each block by its smallest eigenvalue; this is not
a disagreement.

The canonical result also contains the independent direct half-Gram lower
bound. All 512 height segments contribute positively and Arb uses no
mean-value fallback.

## Independent pointwise comparison

At two rational cusp points, the implementation compares:

1. direct evaluation of all six components at `z`, `-z`, `Qz`, and `-Qz`;
2. evaluation of the collapsed integral-lattice Fourier field.

Both complex difference balls contain zero. The largest absolute difference
upper bound is `5.07e-51`.

## Exact plateau: no mass subtraction

The Picard horoball `y>1` is embedded modulo its parabolic subgroup: if
`c != 0`, the height transformation and `|c|>=1` give

\[
y(\gamma P)\le 1/(|c|^2y)<1.
\]

Choose the exact `C2` cutoff

\[
\chi_B(y)=
\begin{cases}
0,&y\le1.01,\\
10t^3-15t^4+6t^5,&1.01<y<1.20,\\
1,&y\ge1.20,
\end{cases}
\quad t=(y-1.01)/0.19.
\]

Start with a normalized subordinate family `sum_j phi_j=1` for every other
chart. Use `W_B` with weight `chi_B` and every other candidate with weight
`(1-chi_B) phi_j`. The weights still sum to one, and the exactly glued
section satisfies `U=W_B` on `Omega_B`.
Since `W_B` is oldspace and rotation odd,

\[
\|P_{Q,-}P_{\rm old}U\|_2
\ge\|W_B\|_{L^2(\Omega_B)}
>0.20380681000849094.
\]

Thus the projected-mass input is admissible with no automorphization
subtraction. What remains is the cancellation-preserving interval residual on
the transition and core: value/first-gradient overlaps and the partition
constants required by Theorem D(K).
