# Current gravity status

Updated: 2026-08-20

This is the authoritative resume bookmark for the active gravity route.  Read
`CLAUDE.md`, this file, `git status --short`, and the commits listed below
before starting another calculation.  Do not infer acceptance from an
uncommitted artifact, a dimension match, or a modular rank alone.

## Last accepted checkpoint

Commit `438dca3` freezes the coordinate-free resolution of the last symbolic
gap in the complete scale--strut carrier.  The resolved artifact SHA-256 is
`ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179`.
The targeted verifier passed `11/11`; no full suite was run.

**DERIVED EXACT; adversarially corroborated; kinematic only.**  On the real
domain `lambda!=1`, `tau!=0`, `(lambda-1)^2-3 tau^2!=0`, the pivot-free wedge
ideal fixes the four disclosed endpoint coefficients exactly.  No additional
real denominator or rank stratum remains.  Both formerly suspect linear
strata pass direct exact rebuilds, and the extra quadratic is the strictly
positive connection norm certificate.

Combined with the three complete finite 600-cell systems, exact rank-240
proof and resolved high-precision audit, this accepts the `1560 x 240`
vertex-scale plus strut map as the complete kinematic boundary-data carrier
in the frozen coordinate choice.  It is not yet a canonical phase space or
a dynamical solution.

Primary files:

- `reproducible/verify_gravity_600cell_corrected_strut_carrier.py`
- `reproducible/verify_gravity_600cell_corrected_strut_alignment.py`
- `reproducible/verify_gravity_600cell_corrected_strut_alignment_adversarial.py`
- `reproducible/verify_gravity_600cell_corrected_strut_canonical_intersection.py`
- `reproducible/verify_gravity_600cell_corrected_strut_canonical_intersection_adversarial.py`
- `docs/gravity/gravity_600cell_corrected_strut_canonical_intersection_result.md`
- `reproducible/verify_gravity_600cell_full_scale_strut_carrier.py`
- `reproducible/verify_gravity_600cell_full_scale_strut_precision.py`
- `reproducible/verify_gravity_600cell_full_scale_strut_symbolic_adversarial.py`
- `reproducible/verify_gravity_600cell_full_scale_strut_symbolic_gap_resolution.py`
- `docs/gravity/gravity_600cell_full_scale_strut_symbolic_gap_resolution_result.md`

## Earlier falsification, now bypassed rather than erased

The tempting interpretation

```text
240 = 120 vertex-scale modes + 120 vertex-lapse modes
```

was tested by exact substitution, not accepted from its dimension.  It fails
on all 3600 face rows for each of the two representatives, with the same exact
first residual under two different local right-inverse graphs.

**DERIVED NEGATIVE.** That particular *local cell-flex lift* of
`Q^120 scale + Q^120 lapse` is not the complete admissible augmented carrier.
The later carrier uses the geometry-derived scale--strut response and does
not revive the refuted lift.

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
- **DERIVED EXACT, adversarially corroborated.** Compatible two-cell geometry
  uniquely fixes the generic scale--strut endpoint response over the stated
  real nondegenerate domain.
- **DERIVED KINEMATIC.** The resulting global `1560 x 240` carrier has exact
  rank 240, correct support and symmetry, and its apparent binary64 Gram
  discrepancy is fully explained by normal-equation conditioning.

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

The complete 240-column scale-plus-strut carrier is now accepted as
kinematics.  The word “dynamic” remains forbidden until an action derivative
is applied.

The immediate calculation is its target-blind pullback through the frozen
Regge action Hessian/strong-equation map.  It must first report ranks,
nullities, symmetries and `2T` sector multiplicities.  Continuum
scalar/vector/tensor labels may be loaded only after that artifact is frozen.
This is the first calculation that can show whether coupled scale-strut data
contain genuine canonical evolution.

The result still does not select a tick, `c`, `G`, a Planck scale or a particle
mass.

## Resume discipline

Run only the verifier for the active mission unless the user explicitly asks
for the full suite.  The latest scoped registry audit has 356 distinct
registered verifiers, no duplicates or unregistered files, plus two
deliberate exclusions; no full suite was run.
Before stopping, update this file with the last accepted commit, the current
open question, and the exact next test.
