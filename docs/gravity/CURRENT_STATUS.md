# Current gravity status

Updated: 2026-08-19

This is the authoritative resume bookmark for the active gravity route.  Read
`CLAUDE.md`, this file, `git status --short`, and the commits listed below
before starting another calculation.  Do not infer acceptance from an
uncommitted artifact, a dimension match, or a modular rank alone.

## Last accepted checkpoint

Commit `6967ed7` freezes a second, mechanically different construction of the
exact 240-direction carrier: one universal rational 6-by-8 response block
repeated over all 600 tetrahedral cells.  Its prior-art gate is `fae9940`,
protocol `dd302d8`, first source `ac5f80f`, recorded first control failure
`0ed3880`, and preregistered coordinate correction `592266a`.  The artifact
SHA-256 is
`0a569e48189c56bc081efcee33f7826fedd52afb93b6135ddb2fec385b56fbdf`.

For both rational representatives, both exact right-inverse graphs, reversed
faces, odd relabelling, and reversed metric sign, the local construction gives:

```text
affine rank                       48 / 48
nonzero local-block entries       48 / 48
nonzero direct residual rows       0 / 6000
one-row-corrupted image           rejected
support of each datum             exact 20-cell vertex star
```

**DERIVED over Q.** The boundary-data carrier is exactly 120 unsigned
vertex-scale edge variations plus 120 arbitrary strut variations.  Every
basis datum has a unique exact cell-flex response satisfying all complete
face equations, and that response is strictly local on vertex stars.

Primary files:

- `reproducible/verify_gravity_600cell_local_data_lift.py`
- `reproducible/gravity_600cell_local_data_lift.json`
- `docs/gravity/gravity_600cell_local_data_lift_protocol.md`
- `docs/gravity/gravity_600cell_local_data_lift_result.md`

## Most recent falsification, now reconciled

The tempting interpretation

```text
240 = 120 vertex-scale modes + 120 vertex-lapse modes
```

was tested by exact substitution, not accepted from its dimension.  It fails
on all 3600 face rows for each of the two representatives, with the same exact
first residual under two different local right-inverse graphs.

**DERIVED NEGATIVE.** That particular *local cell-flex lift* of
`Q^120 scale + Q^120 lapse` is not the complete admissible augmented carrier.
The boundary-data image itself was not refuted by this calculation.

The first execution was preserved in commit `829bf29`.  Two verdict-only
classifier corrections were preregistered and committed separately; no
scientific matrix or residual was changed.  The final targeted run passed
10/10 with the refuted outcome and was frozen in commit `96b912d`.

Primary files:

- `docs/gravity/gravity_600cell_canonical_data_carrier_first_failure.md`
- `docs/gravity/gravity_600cell_canonical_data_carrier_classifier_correction.md`
- `reproducible/verify_gravity_600cell_canonical_data_carrier.py`
- `reproducible/gravity_600cell_canonical_data_carrier.json`

## What the gravity route has established

- **DERIVED.** The carrier is the boundary of the regular 600-cell, a
  triangulated 3-sphere with `f=(120,720,1200,600)`.
- **DERIVED.** The homogeneous cellular weak-lapse Regge equations give the
  all-step discrete Friedmann relation, but do not select an absolute tick.
- **PATTERN with exact finite controls.** Canonical projected refinements move
  the homogeneous acceleration toward the dust-FLRW value `-1/2`.
- **DERIVED.** A schedule-free local flat tetrahedral frustum has six shape
  flexes; variable face transitions are necessary for global closure.
- **DERIVED.** On unequal-scale global slabs the fixed-data system has full
  rank 3600.  On the static slab its 119-dimensional kernel is exactly the old
  vertex-gradient/prism-shift carrier, not a new physical sector.
- **DERIVED NEGATIVE.** The first local lift of the natural 240-dimensional
  vertex scale/lapse data is not compatible globally.
- **DERIVED (modular).** The boundary-data image itself is exactly unsigned
  vertex-edge scale data direct-sum arbitrary strut data over both frozen
  finite fields.
- **DERIVED over Q.** A mechanically different exact solver constructs the
  unique cell-flex lift and verifies zero residual on every face equation.

## What is not yet derived

- **OPEN:** the admissible symplectic/canonical phase space;
- **OPEN:** an action Hessian restricted to that phase space;
- **OPEN:** propagating tensor modes and a dispersion relation;
- **OPEN:** a dynamically selected lapse or physical tick duration;
- **OPEN:** an effective universal speed `c`;
- **OPEN:** Newton's constant `G`, an absolute length, and the Planck scales;
- **OPEN:** any particle mass formula or Standard-Model derivation from this
  gravity calculation.

## Active result and next calculation

The support pattern is now independently proved.  Its prior-art gate is
`fae9940`, protocol `dd302d8`, first source `ac5f80f`, recorded first control
failure `0ed3880`, coordinate correction `592266a`, and first complete
artifact `6967ed7`.  The artifact SHA-256 is
`0a569e48189c56bc081efcee33f7826fedd52afb93b6135ddb2fec385b56fbdf`.

The corrected targeted verifier passed 13/13.  It solved 51,320 global affine
constraints on one universal 6-by-8 block, independently of the prior global
3600-variable elimination.  Every construction gives:

```text
rank                         48 / 48
all block entries nonzero    48 / 48
global residual rows          0 / 6000
corrupted image              rejected
```

**DERIVED.** Every flex row is supported on exactly the four scales and four
struts at its tetrahedron's vertices, and every data column is supported on
exactly its 20-cell vertex star.  The global lift is a repeated exact local
block, not merely a sparse global inverse.

The failed old local formula is also reconciled exactly.  After changing from
its `(sigma,nu)` basis to the accepted `(sigma,raw strut)` basis, it differs
from the selected block by a rank-three local Poincare correction.  Both give
the same local length data, but only the new representative glues globally.
Physical responses agree across the two right-inverse graphs.

The immediate calculation is now the exact Regge boundary action/Hessian
pulled back through this local 240-direction carrier.  It requires a new
prior-art gate and target-disclosed convention protocol before any spectrum
is inspected.  Its job is to decide the constraint/gauge/dynamical split,
especially whether the 120 arbitrary strut data remain lapse freedom or are
selected by the action.

The current result does not select a tick: every strut variation is
kinematically allowed, and no action or constraint equation has yet been
applied to choose it.

## Resume discipline

Run only the verifier for the active mission unless the user explicitly asks
for the full suite.  The last independently reported full baseline was 79/79,
but it is not evidence for later files.  Before stopping, update this file with
the last accepted commit, the current open question, and the exact next test.
