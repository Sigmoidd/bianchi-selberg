# Adding a field or congruence group

The shared engine now separates three contracts: exact group identity,
complete arithmetic inventory, and verified analytic geometry/cusp/scattering
terms. d=2,7,11 are working examples. d=19 uses the same engine once its
arithmetic proofs are supplied. A congruence group uses the same contracts
with its own analytic backend; the level-one formulas cannot be silently
reused for it.

## Public interface

| Object | Responsibility | Location |
|---|---|---|
| `GroupKey` | Exact field, subgroup family, and level ideal | `groups/identity.py` |
| `LevelIdeal` | Canonical integral ideal HNF, norm and exact membership | `groups/identity.py` |
| `GroupData` / `EllipticClass` | Class records, exact witnesses, proof id, analytic backend id | `groups/data.py` |
| `InventoryBackend` / `GroupRegistry` | Bind a self-contained replay to an exact group key and proof id | `groups/registry.py` |
| `QuadraticRing`, `MatrixOps`, `QuadraticOrder`, `Radical` | Reusable exact arithmetic and comparison | `groups/arithmetic.py`, `matrix.py`, `relative_orders.py`, `exact.py` |
| `AnalyticBackend` / `AnalyticRegistry` | Validate scope, prove volume/systole, enclose cusp and scattering terms | `core/backends/` |
| `evaluate`, `certificate_payload` | Common support gate, identity/elliptic sums, assembly and export | `core/assemble.py`, `certificate.py` |
| Uniform command | One CLI for all registered groups | `examples/group_certificate.py` |

Simple full-group callers remain compatible:

```python
from core.assemble import evaluate
from groups import get_group

G = get_group(2)            # exact key: full group, d=2, unit ideal
E = evaluate(G, verbose=False)
S = evaluate(7, include_elliptic=False, verbose=False)  # mechanical screen
```

The old `i`, `omega`, integer d selectors and `bianchi_omega_arb.py`
entry point remain available. The `omega` selector names the field; new
matrix/ideal records always use the **standard** integral basis described
below, whereas the historical flip script explicitly uses its legacy ω
basis.

## Group and level identity

For d=1,2 the standard generator is g=√−d. For d=3,7,11,19,43,67,163
it is g=(1+√−d)/2. An entry (a,b) means a+b g. The corresponding quadratic
relation is available as `QuadraticRing(d).relation`.

`LevelIdeal(d,(a,b,c))` represents columns (a,0), (b,c), with a,c>0 and
0≤b<a. Construction checks closure under multiplication by g, so an
arbitrary rank-two lattice is not accepted as an ideal. The determinant
a c is its norm. `principal` derives the HNF using exact Bezout arithmetic;
`rational(d,n)` gives nO. Associate generators give the same canonical key.

```python
from groups import GroupKey, LevelIdeal

P = LevelIdeal.principal(7, (0, 1))     # (g), norm 2
Q = LevelIdeal.principal(7, (1, -1))    # (1-g), also norm 2, different ideal
A = GroupKey.congruence("gamma0", P)
B = GroupKey.congruence("gamma0", Q)
assert A != B
L = GroupKey.congruence("principal", LevelIdeal.rational(2, 3))
```

The supported families are `full`, `gamma0`, `gamma1`, and `principal`.
For SL lifts, membership means respectively all determinant-one integral
matrices; c∈I; c∈I and a,d≡the same ±1; or the whole matrix congruent to
±I. The sign ensures the definition descends to PSL. Unit-level families
canonicalize to `full`; a nonunit ideal cannot be assigned to `full`.
Groups are never identified by the level norm alone: distinct split prime
ideals often have the same norm.

Creating a key establishes identity and matrix membership. It does not
supply an inventory, cusp data, subgroup index or scattering formula.

## Add d=7,11,19 with one arithmetic adapter

Each field module supplies a `GroupData` factory and a replay function.
Use `groups/builtins.py:d2_group` and `d2_backend` as the concrete pattern.
The d=2-specific archimedean formulas stay in `groups/d2_inventory.py`;
quadratic-order multiplication, norms, discriminants and matrix operations
have moved to shared modules. For example:

```python
from groups.arithmetic import QuadraticRing
from groups.relative_orders import QuadraticOrder

ring = QuadraticRing(7)
trace_one_order = QuadraticOrder(ring, trace=(1, 0), norm=(1, 0))
# Basis (1,g,x,g*x), x²-x+1=0; all operations are exact.
```

The factory sets `inventory_status="self-contained"`, a versioned
`inventory_proof_id`, and `normalization_status="proved"` on its classes
only after the arithmetic proof closes. Each non-cuspidal class supplies
its SL representative, primitive translation, optional endpoint flip,
finite-centralizer size and exact norm triple a+b√c, with optional positive
integer `norm_denominator` (default 1). The primitive norm is then
(a+b√c)/norm_denominator. A cuspidal class can
omit the loxodromic translation/norm; its centralizer and orbital term
must be covered by the cusp backend.

The replay accepts the **actual** GroupData supplied to assembly. It
compares its records against the derived complete inventory and checks
all field-specific facts; returning a separate expected inventory while
ignoring the input is insufficient. Its JSON-compatible dict includes
`proof_id` equal to the backend id. Complete lattice classification,
SL/PSL multiplicity, units, primitive lengths and full orbital normalization
remain mathematical proof obligations; common witness checks do not
replace them. See [INVENTORY_PROOF.md](INVENTORY_PROOF.md).

Once a field module has these proved functions, register it:

```python
from groups import GroupKey, InventoryBackend, default_registry
from groups.builtins import d7_group as make_group
from groups.d7_inventory import verify_group_records  # implemented arithmetic replay

key = GroupKey(7)
backend = InventoryBackend(
    key=key,
    proof_id="d7-arithmetic-v1",
    documentation=("docs/D7_INVENTORY_PROOF.md",),
    replay=verify_group_records,
)
def factory():
    return make_group()  # key/proof id above; trace_backend_id="level-one-v1"

default_registry().register(key, factory, backend, replace=True)
```

`replace=True` explicitly replaces an existing registration (d=7 is already
registered by the builtins; d=11 is also registered; d=19 still has a mechanical-only scaffold).
Ordinary duplicate registration is rejected. Replacing a factory clears
its old proof binding. Keep this initialization in an adapter module or
add it to `register_builtins` once the proof is part of the repo; it needs
no field switch in assembly or certificate export. An isolated
`GroupRegistry` can instead be passed to both `evaluate(...,
group_registry=registry)` and `certificate_payload(...,
group_registry=registry)`.

The common inventory layer checks that the key/proof id match, class labels
are unique, representatives and translations belong to the exact subgroup,
SL determinants are one, elliptic orders hold in PSL, translations commute,
and supplied flips anticommute on the lifts. The registered replay then
proves completeness and binds the arithmetic parameters. Registration is
trusted proof code supplied by the researcher, not ingestion of an
untrusted completeness flag.

## Add a congruence-level analytic backend

A proper subgroup requires its own subgroup index/volume proof, cusp
orbits and lattices, scattering determinant or matrix, and systole or
proved lower bound sufficient for support. Counting elliptics may also
require subgroup conjugacy splitting. Full-group counts, unit norms and
cusp constants do not automatically transfer.

Implement the `AnalyticBackend` protocol:

| Method | Required result |
|---|---|
| `validate(group)` | Reject every group outside the backend's proved scope |
| `geometry(group)` | Replay geometry facts; return `Geometry(group.key, volume_ball, systole_ball)` with finite positive Arb enclosures |
| `terms(group, geometry, k, delta, R, include_elliptic)` | Return finite Arb values for `CE`, `Ch0`, `PARg0`, `PSI`, `PHIINT`; optional diagnostic terms such as `prime` |

`CE` aggregates cuspidal elliptics; `Ch0` aggregates the scattering-at-zero
and parabolic h(0) pieces; `PARg0` aggregates lattice constants; `PSI`
aggregates cusp digamma integrals; `PHIINT` encloses the scattering integral.
If useful, the geometry result can subclass `Geometry` to carry its cusp
or scattering context. The existing `LevelOneGeometry` does this for the
one-cusp constants. `systole_trace=None` is allowed when the new backend
uses a different geometric proof.

Set a versioned `backend_id`, document the derivations in `documentation`,
and register an instance in `core.backends.default_registry()` or an
isolated `AnalyticRegistry`. Set the GroupData's `trace_backend_id` to it.
The full-group backend rejects every proper subgroup even for a mechanical
screen. A geometry result with another group's key is rejected by the core.

The engine reuses test functions, exact spline derivatives, Arb support
comparison, identity term, per-class non-cuspidal elliptic sum, positivity
argument, panel quadrature and certificate export. The bound always uses

\[
 B=I+NCE+CE+Ch0+PARg0+PSI+PHIINT-h(i).
\]

A backend supplies its group-specific analytic blocks and their complete
tails. It must honor `include_elliptic=False` for exploratory screens;
those screens cannot export a spectral certificate. The current generic
sinc family uses integer k≥2 and a finite integration cutoff R≥7.
No congruence-level backend is certified by this refactor: the architecture
and exact membership are ready, and its analytic formulas remain new work.
The existing FEM congruence proofs stay in their separate implementation.

## Commands and certificates

```sh
python examples/group_certificate.py 2 --output /tmp/d2.json
python examples/group_certificate.py 7 --mechanical --output /tmp/d7-screen.json
python -m unittest discover -s tests -v
```

Once registered, a congruence key is selected with `--subgroup gamma0` and
`--level-generator A B`. The CLI deliberately fails if that key has no
factory or analytic backend. Nonprincipal ideal keys can be supplied by
the Python API as explicit HNF.

Schema v2 records `group_identity`, `inventory_backend`, `analytic_backend`,
volume and systole balls, proof sources, exact witness manifest and numerical
parameters. Historical v1 artifacts remain frozen. The new d=2 v2 artifact
has exactly the v1 arithmetic manifest, δ, integration cutoff and endpoint
balls. The tests exercise real d=2 dispatch plus **synthetic, explicitly
noncertifying** congruence fixtures, level membership and backend isolation.

The implemented d=7 adapter is in `groups/builtins.py` and its proof replay
is `groups/d7_inventory.py`. Exact finite quotient and adjacency interfaces
are documented separately in [FINITE_QUOTIENTS.md](FINITE_QUOTIENTS.md).

The d=11 adapter follows the same interface in `groups/builtins.py`. Its
proof and replay are [D11_INVENTORY_PROOF.md](D11_INVENTORY_PROOF.md) and
`groups/d11_inventory.py`; [D11_REVIEW_GUIDE.md](D11_REVIEW_GUIDE.md) records
the finite parameter comparison and regression scope.
