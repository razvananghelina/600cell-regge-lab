# Current gravity status

Updated: 2026-08-19

This is the authoritative resume bookmark for the active gravity route.  Read
`CLAUDE.md`, this file, `git status --short`, and the commits listed below
before starting another calculation.  Do not infer acceptance from an
uncommitted artifact, a dimension match, or a modular rank alone.

## Last accepted checkpoint

Commit `16ae791` freezes a mechanically independent exact-rational lift of
the complete 240-direction canonical boundary-data carrier.  Its prior-art
gate is `4de7e8f`, protocol `1cbe2e1`, and pre-execution registered source
`c02ad2b`.  The artifact SHA-256 is
`1b6ac46a0ea4889f476cc71d51ac464c27caa6d4b6a9b2f2d74ff93da77b123f`.

For both rational representatives, both exact right-inverse graphs, reversed
faces, odd relabelling, and reversed metric sign:

```text
exact pivots of F                 3600 / 3600
candidate consistency             true
nonzero direct residual rows          0 / 6000
one-row-corrupted image rejected   true
exact lift nonzeros               28800
```

**DERIVED over Q.** The boundary-data carrier is exactly 120 unsigned
vertex-scale edge variations plus 120 arbitrary strut variations.  Every
basis datum has a unique exact cell-flex response satisfying all complete
face equations.

Primary files:

- `reproducible/verify_gravity_600cell_rational_data_lift.py`
- `reproducible/gravity_600cell_rational_data_lift.json`
- `docs/gravity/gravity_600cell_rational_data_lift_protocol.md`
- `docs/gravity/gravity_600cell_rational_data_lift_result.md`

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

The exact lift has a rigid support census in all seven constructions:

```text
every flex coordinate uses 8 data coordinates;
every vertex-scale datum touches 20 tetrahedral cells;
every strut datum touches 20 tetrahedral cells.
```

Those are precisely the incidence counts of four scale plus four strut data
per tetrahedron and the 20-cell star of each 600-cell vertex.

**PATTERN.** The current artifact proves the counts but not yet membership of
each support set.  Do not upgrade this to a local formula by numerology.

The immediate calculation is a target-disclosed exact support test: require
each data column's cell support to equal its vertex star and each flex row's
data support to equal the eight coordinates at its cell's four vertices.  If
it passes, extract the exact local 6-by-8 response blocks and reconcile them
with the refuted old local formula.  Only after that should the Regge boundary
action/Hessian be restricted to the accepted carrier.

The current result does not select a tick: every strut variation is
kinematically allowed, and no action or constraint equation has yet been
applied to choose it.

## Resume discipline

Run only the verifier for the active mission unless the user explicitly asks
for the full suite.  The last independently reported full baseline was 79/79,
but it is not evidence for later files.  Before stopping, update this file with
the last accepted commit, the current open question, and the exact next test.
