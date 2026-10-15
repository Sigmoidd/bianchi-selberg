# Picard level-one quasimode: 2026-10-15 verification checkpoint

**No spectral enclosure is claimed.** The certificate now includes:

1. Exact reduction from 14 to **9 active** Ford-cell cosets.
2. Exactly **real, t-reflection-even, Q-odd** 40-orbit projection, correcting half-square Parseval applicability.
3. Independent exact-rational Simpson bounds and two MPFR replays proving `||f||_{1<=y<=5/3} > 1.3e-5` for the projected finite seed.
4. Exact volume bounds `vol(F_core) <= 73/300`, `vol(strip)=4/25`.
5. Reproducible numerical only root-split Gauss and Sobol checks of the reduced analytic residual.

**Two unproved uniform gates** close a spectral enclosure of width `1e-5`: max active automorphy defect `<7e-6` on strip, and sup residual `<2e-10` on Ford core. Neither is certified and no spectral theorem flags were changed.

Read `RESIDUAL_INTERVAL_CLOSURE_GATE.md` first. Replays from this directory:

```sh
python reflection_even_trial.py
python proof_gate_audit.py
python ford_active_cosets.py
python audit_seed_slab_simpson_bounds.py
gcc -O2 -std=c11 -DNT=1024 -DNY=32 -o certify_seed_mass certify_picard_seed_slab_mpfr.c -l:libmpfr.so.6 -lgmp -lm
./certify_seed_mass
# optional independent-resolution check
gcc -O2 -std=c11 -DNT=1536 -DNY=40 -o certify_seed_mass_fine certify_picard_seed_slab_mpfr.c -l:libmpfr.so.6 -lgmp -lm
./certify_seed_mass_fine
OPENBLAS_NUM_THREADS=1 python -m unittest -v test_residual_support_reduction
OPENBLAS_NUM_THREADS=1 python root_split_diagnostics.py
OPENBLAS_NUM_THREADS=1 python sample_uniform_gates.py
```

Python numerical diagnostics require numpy, mpmath, scipy. The directed-rounding proof uses libmpfr6, libgmp and a C compiler. On normal systems with MPFR headers installed, link `-lmpfr`; the source contains an LP64 ABI fallback for this specific environment.
