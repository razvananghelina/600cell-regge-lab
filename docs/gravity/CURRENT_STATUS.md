# Current gravity status

Updated: 2026-08-19

This is the authoritative resume bookmark for the active gravity route.  Read
`CLAUDE.md`, this file, `git status --short`, and the commits listed below
before starting another calculation.  Do not infer acceptance from an
uncommitted artifact, a dimension match, or a modular rank alone.

## Last accepted checkpoint

Commit `f2a35f6` freezes a mechanically different adversarial audit of the
geometry-selected corrected strut carrier against the frozen canonical and
dynamic targets.  Its primary target-blind carrier is `dab941b`, primary
alignment artifact `7ef7a7b`, adversarial protocol `19ea7d3`, and adversarial
verifier registration `2c72165`.  The adversarial artifact SHA-256 is
`3b0fd6da76195279f1beac540c326c61eff5e3172a63bb89baf69502254c5b1f`.

**DERIVED COMPUTATIONAL, adversarially corroborated.** All 14 corrected
parity-sector carrier subspaces differ from the canonical pole-Schur lift.
The 42 fixed comparisons are all separated, with projector distances
`0.997794964..0.998315953`.  A polar/projector/direct-eigenvector audit agrees
with the primary QR/SVD/Schur route to `1.353e-13` and passes 14/14 controls.

The complete extreme-branch interpretation remains **OPEN** because the
last one-dimensional sector has modulus gap `1.006134 < 2`.  The large
projector distance proves non-equality, not zero intersection or global
near-transversality.

Primary files:

- `reproducible/verify_gravity_600cell_corrected_strut_carrier.py`
- `reproducible/verify_gravity_600cell_corrected_strut_alignment.py`
- `reproducible/verify_gravity_600cell_corrected_strut_alignment_adversarial.py`
- `docs/gravity/gravity_600cell_corrected_strut_alignment_result.md`

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

## What is not yet derived

- **OPEN:** the admissible symplectic/canonical phase space;
- **OPEN:** an action Hessian restricted to that phase space;
- **OPEN:** the intersection dimension of corrected geometric struts with
  the canonical strong-equation graph;
- **OPEN:** propagating tensor modes and a dispersion relation;
- **OPEN:** a dynamically selected lapse or physical tick duration;
- **OPEN:** an effective universal speed `c`;
- **OPEN:** Newton's constant `G`, an absolute length, and the Planck scales;
- **OPEN:** any particle mass formula or Standard-Model derivation from this
  gravity calculation.

## Active result and next calculation

The corrected strut formula is target-blind and exact, but its tempting
`119+1` identification with the known dynamic split has failed decisively.
Kinematic face-gluing compatibility and canonical stationarity select
different subspaces.

The immediate calculation is the sectorwise nullity of

```text
G_corrected - C_canonical.
```

Both graphs have literal identity on the same five pole positions, so this
kernel is exactly their intersection coefficient space.  It must be frozen
with rank tolerances and known synthetic controls before execution.  This
will decide how many, if any, corrected pure-strut directions also satisfy
the canonical strong equations.  The complete 240-dimensional
scale-plus-strut action/Hessian pullback follows only after this intersection
census.

The result still does not select a tick, `c`, `G`, a Planck scale or a particle
mass.

## Resume discipline

Run only the verifier for the active mission unless the user explicitly asks
for the full suite.  The latest scoped registry audit has 350 distinct
registered verifiers plus two deliberate exclusions; no full suite was run.
Before stopping, update this file with the last accepted commit, the current
open question, and the exact next test.
