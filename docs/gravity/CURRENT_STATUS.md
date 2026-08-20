# Current gravity status

Updated: 2026-08-20

This is the authoritative resume bookmark for the active gravity route.  Read
`CLAUDE.md`, this file, `git status --short`, and the commits listed below
before starting another calculation.  Do not infer acceptance from an
uncommitted artifact, a dimension match, or a modular rank alone.

## Last accepted checkpoint

Commit `d7fb983` consolidates the exact and adversarial resolution of the last
homogeneous canonical-intersection gap.  The result is zero nonhomogeneous
directions and exactly one homogeneous weak-pole line.  The later transversality
test proves that it is not a tangent of the full pole/lapse equations or a
physical tick.

The immediately subsequent registered transversality test, consolidated in
commit `c5197b2`, passed `6/6`.  Its
artifact SHA-256 is
`d8fd2b0cd71d428d6cef5874b0cd6cf0496f174db13471bdb818a0803d182e0a`.
The full pole equation is transverse to the weak line by a determinant certified
`3.19e13` times outside its error bound.  Therefore the weak line is off-shell
and is not a free tick or gauge direction.  The accepted nonstatic endpoint is
preserved and locally isolated at fixed incoming data.

The final logical composition, consolidated in commit `0e3af17`, now closes
the carrier route.  The corrected
artifact SHA-256 is
`964e993fd9078387eab7064537b5f496d46abfcfd77182671bbc0903ec6e29a4`.
All twelve nonhomogeneous weak intersections are zero, while the only two
homogeneous weak lines are removed by the full pole equation.  Hence all
fourteen full-equation parity/sector intersections are zero.

**DERIVED LOGICAL/COMPUTATIONAL NEGATIVE.** The exact 240-column scale+strut
carrier is not a nonzero full-equation physical tangent space.  This is the kill
boundary for the carrier-intersection selection route, not for the existing
homogeneous roots or unrestricted canonical map.

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
- `reproducible/verify_gravity_600cell_homogeneous_pole_transversality.py`
- `docs/gravity/gravity_600cell_homogeneous_pole_transversality_result.md`
- `reproducible/verify_gravity_600cell_full_equation_carrier_no_go.py`
- `docs/gravity/gravity_600cell_full_equation_carrier_no_go_first_failure.md`
- `docs/gravity/gravity_600cell_full_equation_carrier_no_go_result.md`

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
- **DERIVED COMPUTATIONAL.** The complete pole equation is transverse to that
  weak line.  Hence it is an off-shell lapse response, not a full solution
  tangent, free tick or gauge direction.  The already accepted nonstatic root
  remains locally unique at fixed incoming canonical data.
- **DERIVED LOGICAL/COMPUTATIONAL NEGATIVE.** Combining every sector with the
  full pole equation gives zero complete carrier/canonical intersection in all
  14 parity/sector cells.  This carrier cannot select a nonzero physical
  perturbation tangent.
- **DERIVED COMPUTATIONAL / STRUCTURAL:** independently of that carrier no-go,
  four consecutive homogeneous fixed-mass roots and an action-generated
  1,440-dimensional one-step tangent/Jacobi map already exist.  Their absolute
  time unit and physical constraint quotient remain open.
- **DERIVED SOURCE-LEVEL DISTINCTION:** the exact `lambda=1/2` zero-lapse
  boundary of the repository is a static-product momentum homotopy with
  `Delta L/tau=0`.  It is not the velocity-driven null-strut endpoint in the
  published dust or vacuum-Lambda 600-cell evolutions, where a changing scale
  drives the temporal edge null.  A future nonstatic trajectory may still hit
  that known artifact; this remains open and is a refinement diagnostic.

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

The complete 240-column scale-plus-strut carrier route is closed.  The existing
four-step homogeneous map remains a valid minisuperspace witness, but no
physical anisotropic sector is selected on that fixed carrier.  Do not launch
another fitted subspace or kernel search.

The refined-map feasibility census has now passed `8/8` twice with byte-identical
artifact SHA-256
`ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e`.
The disclosed first infrastructure failure is frozen at `fc92015`, its narrow
control repair was preregistered at `70b13cb`, and the corrected verifier is at
`969bafd`.

**DERIVED:** all 24 colour orders produce distinct internal temporal edge sets
on both `K0=P(sd K_600)` and `K1=P(Esd_2(sd K_600))`.  Their cross-diagonal
intersections are empty and unions contain both diagonals over every spatial
edge.  Thus the spatial geometry plus fixed time orientation does not select a
single simplicial temporal carrier.  The frozen outcome is
`REFINED_MAP_SCHEDULE_ELIMINATION_REQUIRED`.

This does not prove that physical boundary dynamics differs.  The next exact
test is all-24 covariance/equality of the effective boundary quadratic action
on `K0`, after schedule-specific internal variables are eliminated on a common
labelled boundary.  It must be sparse or matrix-free and target-blind.  If the
effective operators differ and no independent geometric rule selects or sums
them, the simplicial refined 600-cell dynamics reaches its kill boundary.  If
they agree, temporal-schedule ambiguity is removed, but a separately derived
coarse/fine phase-space transport is still required before comparing spectra.

The inherited `tau0=0.0102` still supplies, rather than derives, the absolute
time scale.

The source-level comparison in
`gravity_600cell_null_strut_prior_art_reconciliation.md` must be retained in
that refinement mission: the known dynamic null condition is a scale-velocity
condition, not the already closed static homotopy boundary.

The result still does not select a tick, `c`, `G`, a Planck scale or a particle
mass.

## Resume discipline

Run only the verifier for the active mission unless the user explicitly asks
for the full suite.  The post-registration static AST audit finds 364 distinct
registered verifiers, zero duplicate registrations, zero unregistered files,
zero stale registrations and two reasoned deliberate exclusions.  No full
suite was run for the refined feasibility census.
Before stopping, update this file with the last accepted commit, the current
open question, and the exact next test.
