# Eisenstein involution: frozen normalization and exact witness

The historical engine uses one non-cuspidal order-2 class, finite-centralizer
order 2, and N(T0)=7+4 sqrt(3). Its coefficient log(N(T0))/8 remains unchanged.
This note records the issue without changing or re-certifying that input.

In O=Z[omega], omega^2=-1-omega, set

\[
R=\begin{pmatrix}0&-1\\1&0\end{pmatrix},\qquad
X=\begin{pmatrix}\omega&\omega^2\\\omega^2&-\omega\end{pmatrix}.
\]

Integer-coordinate multiplication gives det(R)=det(X)=1,
R^2=X^2=-I, and XR=-RX. Thus the two elements commute in PSL, even
though their SL lifts anticommute. The four distinct classes
{I,R,X,RX} form a Klein-four subgroup of the full PSL centralizer of R.
The portable script verifies all 16 products and all four squares exactly.
It uses no floating point embedding in these computations.

```sh
python scripts/feasibility/flip_check.py
python scripts/feasibility/flip_check.py --bound 2
python /path/to/checkout/scripts/feasibility/flip_check.py --repo /path/to/checkout
```

The theorem's definition of E(R), immediately following Friedman's Theorem
4.1.1 (printed p.42), refers to a maximal finite subgroup of C(R). The witness
proves a subgroup of order 4, **not** maximality. Bounded counts of commuting
or anticommuting matrices are not the order of a maximal finite subgroup:
the full centralizer can also contain loxodromic elements and additional
finite-order elements outside a given finite subgroup.

If one can independently justify |E|=4 while keeping the same class and
primitive norm, the displayed trace-formula coefficient becomes log(N)/16.
The NCE contribution is then **halved**, not doubled. At the historical
parameters it would decrease B by about 0.254723, leaving B near 0.2797.
This is a conditional sensitivity calculation, not a corrected certificate.

The orbital integral with endpoint exchange is now derived in
`D2_INVENTORY_PROOF.md` for the new d=2 records. For the frozen Eisenstein
result, remaining work includes determining the full integral centralizer of the actual class,
prove the maximal finite subgroup and primitive loxodromic norm, and check
element-conjugacy multiplicity. A bounded SL-order search does not settle
these questions. Picard trace-one order-3 lifts cannot have a sign flip,
because conjugation preserves trace while negating changes 1 to -1;
the script also checks that contrast.
