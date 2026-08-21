# Protocol: all-schedule boundary cotangent of the refined on-shell seed

Date: 2026-08-21

Prior-art gate commit: `e7a1545`.

This protocol is frozen before evaluating any refined boundary derivative.

## 1. Frozen inputs

Require exactly:

```text
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
reproducible/gravity_600cell_refined_local_curvature_mass.json
  180010a79177ba16620ebea9847443c57a7a6d2d8a3df71ad6ecb83f454ef091
reproducible/gravity_600cell_refined_local_curvature_mass_adversarial.json
  c59890d12bf929c4677dffed1b932ad8c05ab0ac00980be15ba780e62744c28e
reproducible/gravity_600cell_refined_homogeneous_cotangent_lift.json
  93dd857bff3b406e86d41a8a4b05d6441cb0e3e1c11e4f53d098555b1218924b
reproducible/verify_gravity_600cell_refined_homogeneous_cotangent_lift.py
  154081b12f74ed8597a4b72b37a99219d64c0905829da1566e840fc562b1c20c
reproducible/gravity_600cell_dust_regular_lapse_identity.json
  5079428fade247f730ebc07e5e2eae388b48045cd5201e84afb3186bfc248a51
reproducible/gravity_600cell_projected_rank_edgewise_acceleration_blind.json
  2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2
docs/gravity/gravity_600cell_refined_boundary_cotangent_prior_art.md
  12a32f97129e7fa388261ffd5b16fa4463778ec67af82f676d0365d96300a672
```

Require all upstream accepted outcomes, the exact `15/15` and `16/16`
curvature-mass results, the rank-one/nullity-five cotangent pullback, and the
regular coarse static identity.

## 2. Direct all-schedule action calculation

Parse only function definitions from the frozen stationary-fill source and
rebuild `K0` and all 24 schedule combinatorics.  At 100 decimal digits use the
induced static product coordinates.  Set the evaluator's built-in mass to
zero and add the four independently frozen curvature masses only to the
vertical derivatives.

For every schedule require:

- all six direct cross gradients below `1e-60`;
- all four vertical gradients after curvature-mass addition below `1e-60`;
- accepted Lorentzian branch diagnostics;
- finite gravitational action and all 12 boundary gradients.

This rechecks that the boundary derivatives are taken at an internally
stationary seed rather than an off-shell point.

## 3. Complete boundary covectors

For each schedule and orbit `i` define

```text
P_pre,i  = -dS/dlog(q_old,i),
P_post,i = +dS/dlog(q_new,i).
```

Before any scalar comparison, write all 24 old and new six-vectors.  Require
their componentwise spread below `1e-60` and

```text
P_post=-P_pre
```

within `1e-60`.  A scalar pullback match does not repair schedule-dependent
components.

For the lexicographically first schedule and its reverse, independently
finite-difference all 12 boundary log coordinates at steps `1e-12` and
`5e-13`, Richardson combine the two centered differences and require maximum
relative disagreement with the analytic gradients below `1e-36`.

## 4. Homogeneous pullback and coarse comparison

Require the frozen cotangent artifact to provide orbit-total row

```text
R=(2,2,2,2,2,2).
```

Compute

```text
p_pre,fine  = 2 sum_i P_pre,i,
p_post,fine = 2 sum_i P_post,i.
```

With `K_fine` and `M_fine=K_fine/(8*pi)` from the independently corroborated
curvature artifact, require

```text
p_pre,fine  = -tau0*K_fine/2 = -4*pi*tau0*M_fine,
p_post,fine = +tau0*K_fine/2 = +4*pi*tau0*M_fine
```

within `1e-60` absolute error.

Reconstruct without a stored decimal target

```text
zeta=(pi^2*sqrt(2)/50)^(1/3),
epsilon3=2*pi-5*acos(1/3),
K_coarse=720*zeta*epsilon3,
M_coarse=K_coarse/(8*pi),
p_pre,coarse=-tau0*K_coarse/2.
```

Report:

```text
raw_ratio        = p_pre,fine/p_pre,coarse,
curvature_ratio  = K_fine/K_coarse,
mass_ratio       = M_fine/M_coarse,
normalized_ratio = (p_pre,fine/M_fine)/(p_pre,coarse/M_coarse).
```

Require the first three ratios to agree within `1e-60` and the normalized
ratio to equal one within `1e-60`.  Also require the raw ratio to differ from
one by more than `1e-4`; silently reporting only the normalized equality is a
control failure.

## 5. Falsification controls

1. Add `1e-10` to one boundary component and require the schedule/vector gate
   or pullback gate to fail.
2. Add `1e-6*(1,-1,0,0,0,0)` to the selected old covector.  Its scalar
   pullback must remain exactly unchanged while its vector changes.  This
   proves mechanically that the scalar identity cannot establish component
   selection.
3. Construct the uniform vector with the same pullback and report its distance
   from the action vector.  Do not require that distance to be nonzero before
   seeing the action result; the guaranteed kernel perturbation is the formal
   negative control.

## 6. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_BOUNDARY_COTANGENT_CONTROL_FAILED` for provenance, topology,
   on-shell, branch, finite-difference, corruption or scope failure.
2. `REFINED_BOUNDARY_COTANGENT_SCHEDULE_DEPENDENT` if valid schedules produce
   different complete covectors.
3. `REFINED_BOUNDARY_COTANGENT_IDENTITY_REFUTED` if schedule-independent
   vectors fail the fine curvature/mass pullback identity.
4. `REFINED_BOUNDARY_COTANGENT_SELECTED_EXACT_COARSE` if every gate passes
   and the raw fixed-radius ratio is also one within `1e-60`.
5. `REFINED_BOUNDARY_COTANGENT_SELECTED_RENORMALIZED` if every gate passes,
   normalized momentum agrees exactly, and the raw ratio equals the frozen
   finite-curvature/mass ratio rather than one.

Outcome 5 is **DERIVED COMPUTATIONAL / STRUCTURAL** selection of one refined
action covector, but not exact bare coarse/fine transport at fixed radius.
It points to improved/perfect-action renormalization.  No outcome licenses a
Hessian, mode spectrum, physical tick, `c`, `G`, Planck scale or particle
mass without a separate protocol.

## 7. Deliverables

Register the verifier before execution, run only it twice with byte-identical
JSON and run the static registry audit.  Do not run the full suite or nested
root census.  A significant positive result still needs a mechanically
independent adversarial boundary-gradient audit before consolidation.
