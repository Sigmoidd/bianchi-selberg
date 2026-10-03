# Systoles from traces, with an exhaustive cutoff

Let O be the full ring of integers of Q(sqrt(-d)). Use g=sqrt(-d) for
d=1,2 and g=(1+sqrt(-d))/2 for the odd d in the table. This basis differs
from omega by an integer translation; omega=g-1 when d=3.
Every tau in O is realized as a determinant-one matrix
[[0,-1],[1,tau]]. Conversely every determinant-one matrix has trace in O.
Consequently minimizing loxodromic length over traces gives the systole.

Write tau=x+iy and an eigenvalue as exp((ell+i theta)/2), with ell>0.
Adding its inverse gives
x=2 cosh(ell/2) cos(theta/2) and
y=2 sinh(ell/2) sin(theta/2). Eliminating theta shows

\[
v=\cosh\ell=\frac{A+\sqrt{(A-4)^2+16y^2}}4,\qquad A=x^2+y^2.
\]

Here A and y^2 are rational and are computed with `Fraction`. The larger
root is the loxodromic root; ell=log(v+sqrt(v^2-1)). Real traces with
|tau|<=2 are omitted as elliptic/parabolic. All other listed traces are
loxodromic.

For an eigenvalue of modulus rho>1, |tau|<=rho+rho^-1. Hence |tau|>3
implies ell>2 log((3+sqrt(5))/2), or v>7/2. Each proposed minimum below
has v<=7/2, so every possible shorter trace has |tau|<=3. In the indicated
integral bases this implies |b|<=6 and |a|<=6 for tau=a+b*g. The checker
enumerates that entire box, retains exactly A<=9, identifies equal
algebraic (A,radicand) pairs exactly, and requires a strict Arb inequality
for every other candidate. An inconclusive interval raises an error.

| d | witness tau | cosh(ell0) | ell0 (approximate) |
|---|---|---|---|
| 1 | i | 3/2 | 0.962424 |
| 2 | sqrt(-2) | 2 | 1.316958 |
| 3 | omega | (1+sqrt(21))/4 | 0.862555 |
| 7 | (1+sqrt(-7))/2 | (2+sqrt(32))/4 | 1.265949 |
| 11 | (1+sqrt(-11))/2 | (3+sqrt(45))/4 | 1.534394 |
| 19 | (1+sqrt(-19))/2 | (5+sqrt(77))/4 | 1.907926 |
| 43,67,163 | 3 | 7/2 | 1.924847 |

This replaces the old Eisenstein decimal norm 2.369205407066905; the exact
expression gives a slightly different value. The six-decimal baseline
enclosures remain consistent.

```sh
python -m groups.systoles
```

The same module checks class number one by enumerating primitive reduced
positive definite binary forms (a,b,c) of discriminant D. Reduction gives
|b|<=a<=c, hence |D|=4ac-b^2>=3a^2. Enumerating a<=sqrt(|D|/3), with the
usual nonnegative-b convention at the boundaries, finds exactly one form
for every supported D. This uses the standard form/ideal-class
correspondence and reduction theorem; it is not a classification of
elliptic embeddings or their SL2-conjugacy classes.
