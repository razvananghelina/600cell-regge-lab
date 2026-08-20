# Current gravity status

Updated: 2026-08-20

This is the authoritative resume bookmark for the active gravity route.  Read
`CLAUDE.md`, this file, `git status --short`, and the commits listed below
before starting another calculation.  Do not infer acceptance from an
uncommitted artifact, a dimension match, or a modular rank alone.

## Last accepted checkpoint

Commit `d7fb983` consolidates the exact and adversarial resolution of the last
homogeneous canonical-intersection gap.  The result is zero nonhomogeneous
directions and exactly one homogeneous weak-pole line.  It is not yet a tangent
of the full pole/lapse equations or a physical tick.

The earlier carrier foundation in commit `438dca3` freezes the coordinate-free
resolution of the last symbolic gap in the complete scale--strut carrier.  The
resolved artifact SHA-256 is
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

The primary multiprecision resolver is frozen in commit `99efe9b`, with artifact
SHA-256
`75351ae4dfde26dd75ed8faa927b0a49cd725d83c7629d4545268030b54e2706`.
It passed `17/17`. The registered direct-minor adversarial audit passed `7/7`,
certifying all `48/48` cross-precision nonhomogeneous minors; its artifact
SHA-256 is
`ecf02fd76b0c1d4d95cd206c639a027400c2053bdb1850018d57ff2721861db3`.

**DERIVED COMPUTATIONAL; adversarially replicated.** The complete canonical
scale--strut carrier/action intersection is zero in every nonhomogeneous sector.

The homogeneous gap is now closed.  The primary exact resolver passed `10/10`
and is frozen in commit `3ee5c55`, artifact SHA-256
`70d7583756acdbee77893f98d57054ab074d9353a86247840cc1eb2c7b6be931`.
It derives the generator `sigma=-lambda*p_z, c=p_s` and combines it with all
frozen rank minors to prove D/K nullity one.

The first adversarial run remains an honest `6/7 CONTROL_FAILED` artifact in
commit `5d43620`: its line/rank tests passed but two absolute corruption
thresholds were mis-scaled.  After the disclosed repair preregistered in
`fa798f5`, the fresh P200G direct-matrix reconstruction passed `7/7`, with
artifact SHA-256
`fab74a26ae940cf0e65f26a4f6f167285cc269e282c40d7a630f37d65ba7ab07`.

**DERIVED COMPUTATIONAL; adversarially replicated after disclosed control
repair.** The homogeneous intersection is exactly one-dimensional.  This is a
weak-pole canonical response, not yet a full-equation solution tangent or a
physical tick.

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
- `reproducible/verify_gravity_600cell_full_scale_strut_canonical_intersection.py`
- `docs/gravity/gravity_600cell_full_scale_strut_canonical_intersection_first_result.md`
- `reproducible/verify_gravity_600cell_full_scale_strut_canonical_precision.py`
- `reproducible/verify_gravity_600cell_full_scale_strut_canonical_precision_adversarial.py`
- `docs/gravity/gravity_600cell_full_scale_strut_canonical_precision_primary_result.md`
- `docs/gravity/gravity_600cell_full_scale_strut_canonical_precision_adversarial_result.md`
- `reproducible/verify_gravity_600cell_full_scale_strut_homogeneous_resolution.py`
- `reproducible/verify_gravity_600cell_full_scale_strut_homogeneous_resolution_adversarial.py`
- `docs/gravity/gravity_600cell_full_scale_strut_homogeneous_resolution_result.md`
- `docs/gravity/gravity_600cell_full_scale_strut_homogeneous_resolution_adversarial_first_result.md`
- `docs/gravity/gravity_600cell_full_scale_strut_homogeneous_resolution_adversarial_result.md`

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
- **DERIVED COMPUTATIONAL, adversarially replicated.** All nonhomogeneous
  complete-carrier/action intersections are zero. The primary interval-Gram
  certificates and 48 mechanically different direct-minor certificates agree.
- **DERIVED COMPUTATIONAL, adversarially replicated after disclosed control
  repair.** The homogeneous sector has exactly one weak-pole canonical line.
  Its exact generator, direct D/K matrix kernels, 50/50 rank certificates and
  both parity representatives agree.  The failed first adversarial control is
  retained rather than erased.

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

The complete 240-column scale-plus-strut carrier is accepted as kinematics.  The
action pullback has excluded every nonhomogeneous intersection and selects one
and only one homogeneous weak-pole line by two mechanically different methods.

The exact next calculation must add the pole/lapse equation omitted from the
weak canonical graph.  Preregister the full-equation differential on the now
frozen generator.  If it is nonzero, the surviving weak line is off-shell and
the present carrier has no dynamical tangent.  If it vanishes, construct the
full homogeneous solution tangent and test finite continuation.  Do not infer
gauge or time from the weak line alone.

The result still does not select a tick, `c`, `G`, a Planck scale or a particle
mass.

## Resume discipline

Run only the verifier for the active mission unless the user explicitly asks
for the full suite.  The latest static registry audit has 361 distinct
registered verifiers, no duplicates or unregistered files, plus two deliberate
exclusions; no full suite was run.
Before stopping, update this file with the last accepted commit, the current
open question, and the exact next test.
