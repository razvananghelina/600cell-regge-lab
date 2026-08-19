# Current gravity status

Updated: 2026-08-19

This is the authoritative resume bookmark for the active gravity route.  Read
`CLAUDE.md`, this file, `git status --short`, and the commits listed below
before starting another calculation.  Do not infer acceptance from an
uncommitted artifact, a dimension match, or a modular rank alone.

## Last accepted checkpoint

Commit `6236103` freezes a mechanically different adversarial certificate
that the corrected pure-strut graph has zero intersection with the canonical
strong-equation graph.  The primary intersection artifact is `b64fd83`,
adversarial protocol `2afb0c7`, first serialization failure `1438df3`, and
serialization-only repair `e6fcb66`.  The final adversarial artifact SHA-256
is `c186260ee9520eac59658e3290fb1f4502fd9a7d92f533e8774506cd30e9d03b`.

**DERIVED COMPUTATIONAL, adversarially corroborated.** Every one of the 14
minimal parity-sector differences `G_corrected-C_canonical` has full column
rank.  Hence the global 120-column graphs intersect only at zero for each
staircase parity.  The primary smallest-singular margin is `2.792e6` times
its uncertainty; the independent pivoted-QR/Frobenius lower-bound margin is
`8.536e5`.

Pure struts are kinematically admissible but are not independent canonical
lapse freedom on this curved slab.  This does not remove struts from
evolution: any canonical direction in the complete carrier must mix them
with scale/form data.

Primary files:

- `reproducible/verify_gravity_600cell_corrected_strut_carrier.py`
- `reproducible/verify_gravity_600cell_corrected_strut_alignment.py`
- `reproducible/verify_gravity_600cell_corrected_strut_alignment_adversarial.py`
- `reproducible/verify_gravity_600cell_corrected_strut_canonical_intersection.py`
- `reproducible/verify_gravity_600cell_corrected_strut_canonical_intersection_adversarial.py`
- `docs/gravity/gravity_600cell_corrected_strut_canonical_intersection_result.md`

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
- **DERIVED.** The non-static trapezoid geometry selects a unique rank-120
  pure-strut response with both endpoints on every staircase diagonal.
- **DERIVED COMPUTATIONAL, adversarially corroborated.** That corrected
  carrier is not the full canonical weak lift and is not either frozen
  fixed-count hyperbolic candidate.
- **DERIVED COMPUTATIONAL, adversarially corroborated.** The corrected
  pure-strut carrier and canonical graph have zero intersection globally;
  struts cannot form a nonzero canonical direction without scale/form data.

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

The corrected strut formula is target-blind and exact, but its tempting
`119+1` identification with the known dynamic split has failed decisively.
The stronger intersection census is now complete: no pure-strut direction
survives canonical stationarity.

The immediate calculation is the complete 240-column dynamic
scale-plus-strut carrier pulled through the frozen action Hessian/strong-
equation map.  It must first report target-blind ranks, nullities and `2T`
sector multiplicities.  Continuum scalar/vector/tensor labels may be loaded
only after that artifact is frozen.  This is the first calculation that can
show whether coupled scale-strut data contain genuine canonical evolution.

The result still does not select a tick, `c`, `G`, a Planck scale or a particle
mass.

## Resume discipline

Run only the verifier for the active mission unless the user explicitly asks
for the full suite.  The latest scoped registry audit has 352 distinct
registered verifiers plus two deliberate exclusions; no full suite was run.
Before stopping, update this file with the last accepted commit, the current
open question, and the exact next test.
