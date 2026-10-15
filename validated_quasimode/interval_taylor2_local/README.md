# Picard narrow quasimode: local interval Taylor2 checkpoint

**Do not claim a global spectral enclosure from this archive.**

Source files:
- `seed_interval_prototype.py`: baseline interval boxes; cancellation is not preserved across finite cells.
- `centered_defect_iv.py`: first-order mean-value defect prototype.
- `taylor2_defect_iv.py`: second-order rigorous interval jets for five nonidentity cosets on the triangular wedge.
- `certify_local_mass_tiles.py`: box-local interval mass-transfer correction.
- `local_residual_iv.py`: box-local analytic commutator residual using the defect jets.
- `independent_mp_point.py`: independent direct Bessel value/gradient numerical baseline.
- `test_local_interval.py`: exact 160-to-21 Fourier identity, positive/negative checks, mutation rejection.
- `inputs/`: the upstream frozen projected rational trial, group recurrence, and 24-norm MPFR-validated radial coefficient table and provenance. For upstream verification, also consult PR #14 in `Sigmoidd/bianchi-selberg`.
- `TAYLOR2_LOCAL_INTERVAL_REVIEW.md`: derivation, results, limitations and implementation next step.

Run: `OPENBLAS_NUM_THREADS=1 python -m unittest -v test_local_interval`.
