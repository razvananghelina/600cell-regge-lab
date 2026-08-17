# Nonlinear boundary covariance: quadratic no-go on the frozen quotient

Date: 2026-08-17

Status: post-result audit of the preregistered 32-case nonlinear comparison.
Only the targeted verifier was run; the full suite was not run.

## Provenance

- prior-art gate: `526a202`
- preregistered protocol: `05f76c3`
- frozen Stage-A cases: `b6370bd`
- Stage-B verifier registered before execution: `54c894f`
- frozen input SHA-256:
  `2104c69ba6b21d3a3d92c7071d7f2702cb7d33f7f0e3ff17954f64c469f0c01d`
- result artifact:
  `reproducible/gravity_600cell_dust_nonlinear_boundary_covariance.json`
- result SHA-256:
  `a1e00071fa41f986dfaee84ea6e7689a14c50823f6c87d76889e6cb9346a7e3f`

No continuum answer, speed, full-carrier result or desired nonlinear outcome
was parsed.  The 32 inputs and both linear seeds per input were committed
before any perturbed nonlinear equation was evaluated.

## Complete hypotheses

The result concerns only the following construction.

1. The carrier is the order-24 orbit quotient of the accepted Lorentzian
   600-cell dust slab: 30 old boundary edges, 35 internal edge orbits and 30
   new boundary edges.
2. The two compared slabs are the previously frozen even and odd vertex
   schedules.  They are not isomorphic as ordered four-dimensional slabs.
3. Their boundary coordinates are compared using the unique preregistered
   map tagged `IDENTICAL_PHYSICAL_EDGE_SETS`.
4. The action, dust mass, background and branch conventions are exactly those
   in the frozen audited sources.  No perfect-action or counterterm correction
   is added.
5. The inputs are the four older, target-blind Helmert shape directions, pure
   position and pure pre-momentum rays, both signs and the two frozen levels
   one half and one.
6. Each parity and case is solved twice with the independently calibrated
   operational and validation fixed base Jacobians.  The complete 65 real
   canonical equations are solved; internal variables are not held fixed.
7. The compared output is the dimensionless 60-vector
   `(log q_new, p_post/p_star)`.

No claim below covers the other 25 shape directions, mixed rays, the full
720-edge carrier, matter perturbations, a refinement limit or arbitrary
triangulations.

## Mechanical result

The targeted run passed `8/8` controls and returned

```text
classification counts = BROKEN: 32
outcome = NONLINEAR_BOUNDARY_COVARIANCE_BROKEN_ON_FROZEN_CASES
```

All 128 nonlinear solves converged.  Every action evaluation remained on the
same certified Lorentzian branch.

| direction | sector | half-level defect range | full-level defect range | observed order range |
|---:|---|---:|---:|---:|
| 1 | position | `9.8506040842e-13 .. 9.8506040857e-13` | `3.9402416432e-12 .. 3.9402416444e-12` | `2.0000000035 .. 2.0000000037` |
| 1 | momentum | `1.1583448408e-12 .. 1.1583448411e-12` | `4.6333793732e-12 .. 4.6333793753e-12` | `2.0000000031 .. 2.0000000034` |
| 2 | position | `1.8022366739e-23 .. 1.8022734623e-23` | `7.2088731222e-23 .. 7.2091674292e-23` | `1.9999852760 .. 2.0000147249` |
| 2 | momentum | `2.1093274302e-23 .. 2.1093948123e-23` | `8.4371749579e-23 .. 8.4377140151e-23` | `1.9999769567 .. 2.0000230427` |
| 3 | position | `1.8596808724e-23 .. 1.8597077149e-23` | `7.4386698065e-23 .. 7.4388845462e-23` | `1.9999895884 .. 2.0000104120` |
| 3 | momentum | `2.1765893375e-23 .. 2.1766385391e-23` | `8.7062589482e-23 .. 8.7066525613e-23` | `1.9999836942 .. 2.0000163058` |
| 4 | position | `1.0323664645e-24 .. 1.0323780076e-24` | `4.1294427724e-24 .. 4.1295351170e-24` | `1.9999919347 .. 2.0000080656` |
| 4 | momentum | `1.2042622911e-24 .. 1.2042833486e-24` | `4.8170070506e-24 .. 4.8171755101e-24` | `1.9999873869 .. 2.0000126133` |

All 16 preregistered half/full diagnostics are `QUADRATIC_COMPATIBLE`.
The complete observed-order interval is

```text
1.9999769567352574 .. 2.000023042677222
```

Thus the earlier tangent equality is real, but it is only a first-order
statement.  The first resolved schedule dependence is quadratic in every
frozen ray.

## Hostile numerical audit

- All 128 solves report `CONVERGED` in 12--15 iterations.
- Final residual infinity norms lie in
  `1.6686e-60 .. 3.1253e-56`.
- Final fixed-J correction norms lie in
  `1.0349e-55 .. 2.6843e-52`.
- Output changes under the last correction lie in
  `4.8681e-54 .. 1.4084e-50`.
- The smallest defect-to-uncertainty ratio is `2.3535e26`; the largest is
  `7.6566e40`.  Therefore the classification does not depend on the factor
  `100` in the preregistered `BROKEN` gate.
- The even and odd old boundary, new boundary and pre-momentum backgrounds
  agree exactly at stored precision under the physical-edge permutation.
  Every paired input ray also agrees exactly.  The effect is not a mismatched
  background or input.
- Across all cases the boundary-position part of the defect ranges from
  `1.4584e-31` to `6.5457e-19`, whereas the post-momentum part ranges from
  `1.0324e-24` to `4.6334e-12`.  The resolved breaking is therefore dominated
  by the outgoing momentum.
- Directions 2--4 are suppressed by many orders of magnitude relative to
  direction 1.  They are numerically resolved, but this finite census does
  not establish an exact selection rule or its cause.

The operational/validation difference is a convergence diagnostic, not a
rigorous interval-arithmetic bound.  Nevertheless, the 100-decimal action
evaluation, exact stored-coordinate matching, independent fixed-J solves,
tiny final corrections, enormous defect margins and coherent factor-four
half/full response jointly rule out ordinary floating-point noise as an
explanation.  A formal interval proof remains **OPEN**.

## Post-result primary-source audit

The broad phenomenon is already known.

- Dittrich and Hoehn show that linearized Regge calculus can possess exact
  symmetries while higher-order dynamics breaks them and produces quadratic
  pseudo-constraints: [From covariant to canonical formulations of discrete
  gravity](https://arxiv.org/abs/0912.1817).
- Bahr and Dittrich find that curved Regge solutions generally lack exact
  discrete gauge symmetries: [(Broken) Gauge Symmetries and Constraints in
  Regge Calculus](https://arxiv.org/abs/0905.1670).
- Restoring exact continuum dynamics requires an improved/perfect action in
  the framework of [Improved and Perfect Actions in Discrete
  Gravity](https://arxiv.org/abs/0907.4323).
- Four-dimensional discretization independence is generically obstructed and
  may require nonlocal structure: [Discretization independence implies
  non-locality in 4D discrete quantum gravity](https://arxiv.org/abs/1404.5288).
- Canonical evolution generated by the discrete action is standard in
  [Canonical simplicial gravity](https://arxiv.org/abs/1108.1974).

The search also found established 600-cell dust evolution work, including
[A Parallelizable Implicit Evolution Scheme for Regge
Calculus](https://arxiv.org/abs/gr-qc/9411008) and [The Friedmann universe of
dust by Regge Calculus: study of its ending
point](https://arxiv.org/abs/gr-qc/0009093).  It did not locate this exact
order-24 two-schedule boundary-map census.  A search cannot prove novelty, so
external novelty is **OPEN**.

## Status ledger

- **DERIVED COMPUTATIONAL:** on all 32 frozen rays, the two schedules give
  different nonlinear canonical boundary outputs.
- **DERIVED COMPUTATIONAL:** the difference is resolved far above the
  calibrated numerical uncertainty and is quadratic-compatible in all 16
  half/full comparisons.
- **DERIVED COMPUTATIONAL:** the defect is dominated by outgoing momentum.
- **STRUCTURAL:** the exact tangent covariance is the common first derivative
  of two distinct nonlinear boundary maps, not equality of the maps.
- **STRUCTURAL NEGATIVE:** the bare order-24 Regge-dust action does not define
  schedule-independent nonlinear evolution on the tested neighbourhood.
- **PATTERN:** directions 2--4 show very strong additional suppression; four
  directions are too few to infer a representation-theoretic rule.
- **OPEN:** whether a geometry-derived improved/perfect action restores the
  covariance, whether the defect decreases under genuine refinement, and
  whether the full 720-edge theory behaves similarly.
- **OPEN:** external novelty of the exact finite computation.

## Consequence for the programme

This closes the tempting interpretation of the linear result as nonlinear
triangulation independence.  It is an honest no-go for the bare quotient,
not a no-go for Regge evolution itself.

The present quotient therefore cannot yet supply a unique physical tick:
different legitimate internal schedules agree through first order but predict
different outgoing momenta at second order.  A continuation must either

1. derive a preferred schedule from additional geometry, or
2. derive an improved/perfect effective action and test it without fitted
   counterterms, or
3. demonstrate convergence of the schedule defect under a preregistered
   refinement family.

Until one of these succeeds, claims of a derived causal speed, Planck time or
nonlinear gravitational dynamics from this tick are **OPEN**.
