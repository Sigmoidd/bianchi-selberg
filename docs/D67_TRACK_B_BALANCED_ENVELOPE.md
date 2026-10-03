# Local balance with an exact inherited scalar budget

Use the notation of D67_TRACK_B_VERTICAL_ENVELOPE.md, and set
t=a_min R²/b_min. If t=0 take a_T=a_min and b_T=b_min. Otherwise let p be the
smaller of 1/2 and the upward dyadic rounding of sqrt(t) to denominator 2²⁰.
Define

\[
a_T=\frac{a_{\min}p}{p+t},\qquad b_T=(1-p)b_{\min}.
\]

The reverse-Young proof with reserve p b_min gives this lower enclosure for
every positive p, including the capped case p=1/2. The zero-radius case follows
without Young's inequality.

Compare with a_old=a_min/(1+2t), b_old=b_min/2. Always b_T >= b_old.
If p=1/2 the enclosures are equal. Otherwise p >= sqrt(t), so

\[
\frac{a_{old}}{a_T}=
\frac{1+t/p}{1+2t}\le\frac{1+\sqrt t}{1+2t}\le\frac98.
\]

The last inequality is exact: with x=sqrt(t),
9(1+2x²)-8(1+x)=18x²-8x+1=18(x-2/9)²+1/9 > 0.
Both lower matrices use the same h0. Their inverse metric formulas therefore
give B_T⁻¹ <= (9/8)B_old⁻¹ in quadratic-form order. Every metric diameter,
and hence gamma², inherits this comparison.

Combining with the centroid constant gives

\[
\gamma_{balanced}^2\le
\frac98\frac{4983}{8000}\gamma_{original}^2.
\]

balanced_envelope.py computes this fraction from the frozen original ledger.
It is below 1/11, so the direct-index scalar condition with eta=1/10 passes.
The scalar comparison requires the same exact coefficient intervals and h0
as the original ledger. A diagnostic box/hull relaxation must not silently
replace those intervals in a final certificate without replaying its own
scalar bound.

The finite matrix still needs independent rigorous enclosures and a verified
negative-index bound. This note does not certify spectral exclusion.


A second valid choice is eta=1/13. The same exact inherited gamma bound is
less than 1/14, giving c0=1-14 gamma²>0.0189727. This reduces finite mass
inflation to 14/13, with a smaller positive interpolation-error margin.
The exact budget is recorded separately in balanced_eta13_budget.json;
it retains the same coefficient-interval and h0 requirements.
