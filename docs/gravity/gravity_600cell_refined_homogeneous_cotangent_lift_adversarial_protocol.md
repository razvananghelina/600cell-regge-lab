# Adversarial protocol: homogeneous refined cotangent-lift freedom

Date: 2026-08-21

This Rule-4 audit is frozen after the primary artifact and before constructing
an adversarial witness.

## 1. Frozen inputs

Hash and use exactly:

```text
reproducible/gravity_600cell_refined_homogeneous_cotangent_lift.json
  93dd857bff3b406e86d41a8a4b05d6441cb0e3e1c11e4f53d098555b1218924b
docs/gravity/gravity_600cell_refined_homogeneous_cotangent_lift_first_result.md
  8bc42b3174eebd192dfa97898e61f4f0ea3eef0b91d30848f67183a341907ec4
reproducible/verify_gravity_600cell_refined_homogeneous_cotangent_lift.py
  154081b12f74ed8597a4b72b37a99219d64c0905829da1566e840fc562b1c20c
reproducible/gravity_600cell_refined_canonical_map_feasibility.json
  ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e
```

Require the primary `12/12` outcome, but do not import or execute any primary
function.

## 2. Mechanically different construction

Do not call a matrix-rank, nullspace, SVD or pseudoinverse routine.

For unit coarse momentum, construct six orbit-total lifts

```text
L_i = e_i/2,  i=0,...,5.
```

Every lift must satisfy `2 sum_j L_i,j=1`.  Relative to `L_5`, take the five
differences `d_i=L_i-L_5`.  Their first five coordinates form the diagonal
matrix `diag(1/2,...,1/2)`, whose determinant must be computed directly as
`1/32`.  This proves five independent homogeneous invisible directions
without a rank routine.

For per-edge momenta with populations `N_i`, independently construct

```text
mu_i = e_i/(2 N_i).
```

Their corresponding five-coordinate determinant is

```text
product_{i=0}^4 1/(2 N_i),
```

which must be a nonzero exact rational.

## 3. Exact satisfiability attacks

Use Z3 rational arithmetic, not floating point.

1. Two six-component lifts subject only to the canonical equation must admit
   a model in which they differ.
2. Two one-orbit lifts subject to the canonical equation must be forced equal;
   adding inequality must be `UNSAT`.
3. Impose the additional equations `alpha_i=alpha_j` for all six actual
   orbits.  The unique unit-momentum lift must be `(1/12,...,1/12)`, and two
   distinct such lifts must be `UNSAT`.

The third test is an intentionally stronger control, not an allowed repair.
The frozen refined artifact has four distinct barycentric rank-class sizes
`(120,720,1200,600)`.  Spatial `H4` preserves face dimension, so it fixes the
four rank colours and consequently all six rank-pair orbit labels.  The full
`S6` permutation symmetry needed by the control is absent from the declared
geometry.

## 4. Convention and falsification controls

- Repeat the explicit-lift equations at coarse momentum `-7/3`; the affine
  differences must remain nonzero.
- Reverse the orbit order and repeat the direct determinant tests.
- A deliberately corrupted proposed lift `e_0/3` must fail the unit-momentum
  canonical equation.
- Replacing one positive population by zero must make the corresponding
  per-edge unit lift undefined; no division-by-zero case may be classified as
  convention equivalent.

## 5. Outcomes

Use the first applicable label:

1. `ADVERSARIAL_REFINED_COTANGENT_CONTROL_FAILED` if provenance or a control
   fails;
2. `ADVERSARIAL_REFINED_COTANGENT_UNDERDETERMINATION_CORROBORATED` if all
   explicit witnesses, exact determinants and Z3 outcomes pass;
3. `ADVERSARIAL_REFINED_COTANGENT_DISAGREEMENT_OPEN` otherwise.

Only outcome 2 permits consolidation of the narrow **STRUCTURAL** result.
It still does not forbid an action-selected or perfect-action lift.

Create a separate registered verifier and deterministic artifact.  Run it
twice, then run only the static registry audit.  Do not run the full suite,
the nested `H4` search, a slab solver or a spectrum.
