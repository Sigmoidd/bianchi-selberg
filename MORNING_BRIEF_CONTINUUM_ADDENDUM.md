# Morning brief addendum — continuum defect

**Hard map:** unchanged; `rung4_certified=false`.

`continuum_defect_arb.py` now certifies a continuum supremum for the fixed
six-copy finite trial on a rectangular superset of the Humbert core. It uses
160-bit Arb center values, analytic gradient majorants, exact coset pullbacks,
the actual area-1/2 Picard half-torus for normalization, and an Arb final
quotient.

| subdivision | raw continuum upper | normalized upper |
|---:|---:|---:|
| 4 | 1.04046e4 | 1.77113e5 |
| 8 | 2.96556e3 | 5.04820e4 |

The enclosure is finite and decreases under refinement. It remains roughly
nine orders above the target because the first-order derivative majorant sums
absolute modal contributions and loses the cancellations visible at box
centers. This is a rigorous diagnostic enclosure, not an admissible D(K) eta:
the current theorem is single-cusp, while the certified object is the
two-cusp/six-copy induced representation.

An adversarial sub-agent review validated the exact gluing and generator
actions, the centered Lipschitz argument, the Bessel derivative inequality,
domain coverage, fixed-coefficient semantics, and termwise PDE identity. It
found two errors in the first draft: full-torus normalization and binary64
final division. Both were corrected. Normalization now uses exact half-period
Fourier Gram blocks whose smallest eigenvalues are enclosed with Arb's Rump
algorithm, followed by 256-piece interval integration in height. The final
normalized bound remains an Arb ball through endpoint extraction.

See `CONTINUUM_DEFECT_CERTIFICATE.md` and
`continuum_defect_arb_result.json`.
