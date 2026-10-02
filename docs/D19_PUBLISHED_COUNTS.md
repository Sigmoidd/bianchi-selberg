# d=19 published finite-subgroup cross-check

For PSL₂(O₋₁₉), Krämer's original table (27.26), printed p.189 / PDF
page 197 of the [1980 thesis scan](https://theses.hal.science/tel-00628809/document),
was visually read for the d=19 row. Blank entries mean zero (p.188).
The SL groups counted contain ±I, so their PSL quotients give these counts.

| λ₄′ | λ₄* | λ₄ᵀ | μᵀ | μ₂⁻ | λ₆′ | λ₆* | μ₃ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 2 | 0 | 0 | 1 | 2 |

| PSL finite subgroup | Conjugacy classes | Inventory interpretation |
|---|---:|---|
| C₂ | 1 | One involution element class |
| C₃ | 1 | Its inverse generators merge; one order-3 element class |
| V₄, including nonmaximal copies | 2 | Inside the two tetrahedral classes |
| S₃ | 2 | Inverse-conjugating order-3 normalizers exist |
| A₄ | 2 | Does not add new elliptic element orders |

This agrees with the two-class self-contained proof. The cyclic count
λ₄′+λ₄*=λ₆′+λ₆*=1 is a subgroup count, not automatically an element
count; the norm-minus-one A₃ unit supplies the inverse conjugation.

Krämer's [Theorem 6.8](https://arxiv.org/html/1207.6460v7#S6) independently
predicts S₃ and A₄ existence: 19≡1 mod 3 and 19≡3 mod 8. It predicts
no maximal V₄ because 19≡3 mod 4. These congruence tests concern existence;
the numeric multiplicities above were read from the original table.
Schwermer–Vogtmann's [1983 paper](https://warwick.ac.uk/fac/sci/maths/people/staff/karen_vogtmann/research/1983.0058.pdf)
is restricted to Euclidean imaginary quadratic integers, hence does not
supply a d=19 figure. No d=19 count is attributed to it.

The external cross-check does not verify primitive norms or coefficients;
those are proved and replayed in [D19_INVENTORY_PROOF.md](D19_INVENTORY_PROOF.md).
