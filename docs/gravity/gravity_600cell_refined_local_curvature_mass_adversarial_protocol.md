# Adversarial protocol: direct local curvature and dust-free action derivative

Date: 2026-08-21

Primary-result commit: `2881e4c`.

This gate is frozen before writing or executing the adversarial verifier.  It
must not import the primary verifier or recover a gravitational derivative by
subtracting the committed `P1` dust term.

## 1. Frozen primary and source inputs

Require:

```text
reproducible/verify_gravity_600cell_refined_local_curvature_mass.py
  c54f17708a2678b925cfce96fcfc7d6baaeeaf0577bedbf22b5d0435c069fae6
reproducible/gravity_600cell_refined_local_curvature_mass.json
  180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091
docs/gravity/gravity_600cell_refined_local_curvature_mass_primary_result.md
  82b55ed2918b2db4f83123b23cb2e7534a9f813d5705a1e2ee0ac92e481a0f10
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

Require the primary outcome
`REFINED_LOCAL_CURVATURE_MASS_IDENTITY_CONFIRMED_POST_HOC`, `15/15`, 24
schedules, 96 equations and mass-response rank four.

## 2. Actual-incidence curvature path

Rebuild the 120-vertex 600-cell and its 600 tetrahedra from
`commons.build_600cell`.  Construct the barycentric flag complex directly,
including all 2,640 vertices, 14,400 chambers and the set of 17,040 actual
spatial edges.  Do not use the primary six edge-population totals as the
curvature accumulator.

Assign every actual edge its exact rank-pair length and deficit.  Add half of
`l_e epsilon_e` to each endpoint, then group the resulting 2,640 vertex
values by rank.  Require:

- exact f-vector `(2640,17040,28800,14400)`;
- actual rank sizes `(120,720,1200,600)`;
- zero within-rank curvature spread at `1e-70` relative tolerance;
- endpoint sum equal to the direct edge sum within `1e-70` relative;
- the four independently accumulated fractions agree with the primary
  fractions within `1e-68` absolute.

As a known-answer control, apply the same endpoint accumulator to the
unrefined regular 600-cell.  Every one of its 120 vertices must receive the
same fraction `1/120`.  Dropping one refined endpoint contribution while
retaining the edge total must fail conservation.

## 3. Direct dust-free Lorentzian action path

Parse only function definitions from the frozen stationary-fill source; do
not execute its module body.  Rebuild schedule combinatorics for order
`(0,1,2,3)` and its reverse `(3,2,1,0)`.  Evaluate its accepted Lorentzian
action with `geometry["mass"]=0`, so the returned rank-lapse gradients are
direct gravitational derivatives rather than `P1`-subtracted values.

At 100 decimal digits require for both schedules and all four ranks

```text
dS_grav/dlog(rho_r)=tau0*K_r/2
```

within `1e-60` absolute error.  Add the independently accumulated masses
`K_r/(8*pi)` analytically and require all four total lapse gradients below
`1e-60`.  Adding `P1` masses `M/4` must leave a norm above `1e-4`.

Independently finite-difference the dust-free complete action in each of the
four `rho_r` log coordinates at steps `1e-12` and `5e-13`, Richardson combine
the two centered differences, and require relative agreement with the direct
analytic derivative below `1e-36`.

## 4. Scope and outcomes

The adversarial verifier must record that it imports no primary function,
runs no stationary root search or nested census, and computes no Hessian,
spectrum, continuum target or physical constant.

Use the first applicable outcome:

1. `ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CONTROL_FAILED` for provenance,
   topology, known-answer, negative-control, finite-difference, branch or
   scope failure;
2. `ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_DISAGREEMENT` if both paths are
   valid but either disagrees with the primary identity;
3. `ADVERSARIAL_REFINED_LOCAL_CURVATURE_MASS_CORROBORATED` otherwise.

Only outcome 3 permits the consolidated label **DERIVED COMPUTATIONAL /
STRUCTURAL, adversarially corroborated, post-hoc discovery**.  It still does
not turn the rank density into homogeneous dust or derive a physical tick,
`c`, `G`, Planck units or particle masses.

## 5. Deliverables

Register the verifier before execution, run it twice with a byte-identical
artifact and perform only the static registry audit.  Do not run the full
suite.
