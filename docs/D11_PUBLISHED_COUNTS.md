# d=11 cross-check against published finite-subgroup counts

Scope: Γ=PSL₂(O₋₁₁), full group, discriminant D=−11. This external
cross-check supplements, and does not replace, the self-contained inventory
proof. **Subgroup conjugacy classes are not element conjugacy classes.**

| Nontrivial finite subgroup type | Γ-conjugacy classes | Relation to the inventory |
|---|---:|---|
| C₂ | 1 | Exactly one nonidentity generator, hence one involution element class |
| C₃ | 1 | Two inverse generators, not Γ-conjugate; hence two order-3 element classes |
| V₄ (D₂), including nonmaximal copies | 2 | Each is the normal V₄ in one of the two A₄ classes; these are not additional involution element classes |
| S₃ (D₃) | 0 | No conjugation exchanging the two C₃ generators |
| A₄ | 2 | Maximal noncyclic finite subgroups; not extra elliptic element types |

The resulting **element** inventory is therefore one involution and two
inverse order-3 classes, in agreement with `d11-arithmetic-v1`. No count
mismatch was found. This cross-check does not independently verify the
primitive norms or the orbital coefficient normalization.

## Original Krämer table row

Table (27.26) counts finite subgroups of SL₂(O); the entries relevant here
contain ±I and therefore correspond bijectively to their PSL₂(O) quotients.
On printed p.188 Krämer states that zero entries are left blank. The d=11
row on p.189 reads:

| λ₄′ | λ₄* | λ₄ᵀ | μᵀ | μ₂⁻ | λ₆′ | λ₆* | μ₃ |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 1 | 1 | 2 | 0 | 1 | 0 | 0 |

Thus λ₄=λ₄′+λ₄*=1 and λ₆=λ₆′+λ₆*=1; two tetrahedral
classes, no maximal quaternion/V₄ classes and no binary-dihedral/S₃
classes. The total V₄ count includes the two inside tetrahedral groups.
This directly tabulated row agrees with the formula evaluation below.
[Source identification and row transcription](../tests/mutation/d11/published-counts.json)
record the scan hash, page numbers and numeric entries; the full scan is
not copied into the repository.

## Sources read and substitutions made

1. Norbert Krämer, *Imaginärquadratische Einbettung von Ordnungen rationaler
   Quaternionenalgebren, und die nichtzyklischen endlichen Untergruppen der
   Bianchi-Gruppen*, arXiv:1207.6460v7. Read [§6, Theorem 6.8](https://arxiv.org/html/1207.6460v7#S6),
   [§7, Theorem 7.5](https://arxiv.org/html/1207.6460v7#S7),
   and [§8, Corollary 8.6](https://arxiv.org/html/1207.6460v7#S8).
   A PDF was also downloaded and searched with `pdftotext -layout`.
   The [1980 thesis scan](https://theses.hal.science/tel-00628809/document)
   was downloaded from HAL through the local HTTP client after browser
   retrieval failed. Its body is scanned; the original appendix table
   (27.26), printed p.189 / PDF page197, was rendered and visually read.
   The d=11 row is directly confirmed, as transcribed below. The formulas
   were additionally checked in Rahm's published reproduction, below.
2. Alexander D. Rahm, *Accessing the cohomology of discrete groups above their
   virtual cohomological dimension*, Journal of Algebra 404 (2014), 152–175;
   [arXiv:1112.4262](https://arxiv.org/pdf/1112.4262), §5 (p.16), Appendix A.1
   (pp.19–22), Appendix A.2 (pp.23–26), and references [16], [17], [27].
   It explicitly evaluates Krämer's formulas for subgroup counts. The
   d=11 row of the 3-conjugacy graph table is one circle (p.21), not two
   C₃ subgroup classes. The 2-torsion table includes d=11 in the two-A₄
   component type (p.26). Formula tables, not a recollected count, supply
   the numeric cyclic counts used here.
3. Joachim Schwermer and Karen Vogtmann, *The integral homology of SL₂ and
   PSL₂ of Euclidean imaginary quadratic integers*, Comment. Math. Helv.
   58 (1983), 573–598. Read [author-hosted scan](https://warwick.ac.uk/fac/sci/maths/people/staff/karen_vogtmann/research/1983.0058.pdf),
   §4.1 (pp.580–581) for the finite types and §5.10, figure (p.591) for
   d=11. The diagram shows two vertex orbits with A₄ stabilizers, a C₂
   edge and a C₃ edge, with the stated side identification and no further
   identifications. It agrees with Krämer's two A₄ classes and the cyclic
   counts. It is a cell-stabilizer computation, not a printed table of
   the three elliptic element classes.

The following are explicit **evaluations/inferences from these sources**.
In Krämer's notation t counts prime divisors of D, so t=1. Theorem 6.8
excludes S₃ since 11≡2 mod3 and excludes maximal finite V₄ since 11≡3 mod4.
It admits A₄ since 11≡3 mod8. Theorem 7.5(ii) then gives μ(A₄)=2^t=2.
Corollary 8.6 gives the number of C₂ classes contained in V₄:
2λ₂*=μ(A₄)+3μ(maximal V₄)=2, hence λ₂*=1.

For possible additional C₂ classes, Rahm's Appendix A.2, row m≡3 mod8,
gives λ₄−λ₄*=½(h(Q(√m))−2^(δ−1)). Here δ=1 and h(Q(√11))=1, so that
difference is zero. Thus there is exactly one C₂ class in all of Γ.
The real class number used in this substitution has an elementary check:
the Minkowski bound √44/2=√11<4 restricts ideal norms to ≤3; the prime
above 2 is principal via 3+√11 of norm −2, and 3 is inert. All such ideals
are principal.

For C₃, Rahm's Appendix A.1, row m≡2 mod3, gives λ₆*=0 and
λ₆−λ₆*=(z/2)h(Q(√(3m))). With m=11, the integer 6+√33 has norm 3,
so z=2. Also h(Q(√33))=1: the Minkowski bound √132/2=√33<6 restricts
norms to ≤5; the two primes above 2 are generated by (5±√33)/2 of norm
−2, the prime above 3 by 6+√33 of norm 3, and 5 is inert. Every ideal
in that range is principal. Hence the total C₃ subgroup count is λ₆=1.
These real-field substitutions are separate from the quartic-order class
numbers used by the production inventory proof.

Finally, suppose an element of Γ conjugated a generator R of C₃ to R⁻¹.
In an eigenbasis its determinant-one lift exchanges the two eigenlines,
so its square is −I. In PSL it is an involution, and together with R it
generates S₃. Krämer's S₃ nonoccurrence therefore proves the two generators
cannot be Γ-conjugate. One C₃ subgroup class consequently gives **two**
elliptic element classes. Passing from subgroup counts to element counts
is essential; equating the two counts would undercount the positive term.
