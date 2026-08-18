# Prior-art gate: canonical Legendre rank of the 600-cell dust slab

Date: 2026-08-16

Status: **completed before evaluating the new `65 x 65` canonical-inversion
Jacobian or any of its previously unavailable pre-momentum rows**.

Upstream corrected gluing control: `a766740`

This is a targeted primary-source map, not a novelty proof.  External novelty
remains **OPEN** pending dedicated review.

## 1. Exact new object, carrier and hypotheses

Use one ordered-schedule 600-cell dust slab and the published regular
time-symmetric solution.  The old boundary is held fixed.  In dimensionless
logarithmic coordinates write

```text
y = (log x_1,...,log x_35,
     log q_new,1,...,log q_new,30),
```

where `x` consists of thirty staircase-diagonal squares and five temporal
pole magnitudes, while each boundary coordinate represents one orbit of 24
physical 600-cell edges under the order-24 schedule stabilizer.

Let

```text
e_a = (1/24) partial S / partial log(x_a),            a=1,...,35,
p_pre,i = -(1/24) partial S / partial log(q_old,i),  i=1,...,30.
```

At fixed `q_old` and fixed published `p_target=p_pre`, define

```text
F(y) = (e_1,...,e_35, p_pre,1-p_target,1,...,p_pre,30-p_target,30).
```

The new object is exactly

```text
J_can = partial F / partial y
```

at the published solution, for each of the already selected `even` and `odd`
schedule parities.  It is the derivative of the pre-Legendre inversion
problem

```text
(q_old,p_pre) -> (x,q_new),
```

not the previously computed internal Dirichlet Hessian and not a fitted
evolution operator.

All variables are logarithmic and every residual is a per-edge logarithmic
action derivative, so row and column units are common.  Rank is invariant
under further nonsingular coordinate rescalings, but singular-value
magnitudes are not; the raw coordinate convention must therefore remain
fixed and be printed in the result.

## 2. KNOWN canonical structure

Dittrich and Hoehn formulate discrete actions as Hamilton principal
functions.  Pre/post Legendre maps, their null vectors and resulting
pre/post constraints determine whether a simplicial move defines a regular
canonical transformation or a constrained relation:

- [*Canonical simplicial gravity*](https://arxiv.org/abs/1108.1974).

Their covariant-to-canonical analysis shows that exact vertex-displacement
symmetry on flat Regge backgrounds yields constraints, while curvature can
break the symmetry and produce background-dependent pseudo-constraints:

- [*From covariant to canonical formulations of discrete
  gravity*](https://arxiv.org/abs/0912.1817).

The linearized Pachner-move analysis identifies lapse/shift generators and
propagating lattice-graviton data only after the constraint and gauge
directions have been separated:

- [*Canonical linearized Regge Calculus: counting lattice gravitons with
  Pachner moves*](https://arxiv.org/abs/1411.5672).

An alternative consistent-discretization program explicitly allows the
discrete equations to determine continuum multipliers and replace a
constrained continuum evolution by a regular canonical transformation:

- [*Consistent discretization and canonical classical and quantum Regge
  calculus*](https://arxiv.org/abs/gr-qc/0511096).

Consequently neither an exact lapse null nor a full-rank discrete Legendre
map may be assumed from continuum intuition.  The matrix must decide.

## 3. KNOWN cosmology and convergence boundary

The published dust-filled 600-cell evolution and its five-stage schedule are
not new:

- [De Felice--Fabri 2000](https://arxiv.org/abs/gr-qc/0009093);
- [De Felice--Fabri 2001](https://arxiv.org/abs/gr-qc/0106077).

Collins--Williams-type Regge cosmologies and perturbed lattice universes are
reviewed in [Liu--Williams](https://arxiv.org/abs/1510.05771).  Regge
convergence cannot be inferred from a small equation residual alone; Brewin
and Gentle show why solution convergence and residual convergence can behave
differently and find second-order solution accuracy in a Kasner test:

- [*On the convergence of Regge calculus to general
  relativity*](https://arxiv.org/abs/gr-qc/0006017).

Thus a regular local canonical map would still require refinement and
perturbation tests before being interpreted as continuum gravity.

## 4. Repository controls already established

The following are **CONTROL / DERIVED UPSTREAM** and will not be rediscovered:

1. The complete one-slab action has 35 internal variables, 30 old-boundary
   momenta and 30 final-boundary momenta.
2. On the published Dirichlet problem there is a stationary one-parameter
   collective lapse family.  Its internal tangent is null, while four
   relative lapse modes have nonzero stiffness approximately `4.60497e-8`.
3. After quotienting that internal tangent, the corrected internal Hessian
   has rank 34.
4. Linearized internal stationarity admits 29 zero-sum final-boundary shape
   responses; one compatibility row is homogeneous scale within frozen
   numerical error.
5. The two-slab action and its pre/post sign conventions glue correctly.  Two
   repeated published slabs do not satisfy the shared equation; their cusp is
   nonzero and lies entirely in the homogeneous scale direction.
6. The full boundary contains 720 physical edges.  The present 30 variables
   are orbits of 24 edges, not the unreduced carrier.

These results do **not** determine `rank(J_can)`.  A vector null for the
fixed-two-boundary Dirichlet problem can be lifted when the initial momentum
is fixed, because the pre-momentum may vary along the Dirichlet lapse family.
Calling that vector canonical gauge before testing the bottom 30 rows would
therefore be circular.

## 5. OPEN object and structural prediction

No repository verifier and no located primary source gives the present full
`65 x 65` pre-Legendre Jacobian.  In particular, the mixed block

```text
- partial^2 S /
  [partial log(q_old) partial(log x,log q_new)]
```

has not been evaluated at arbitrary precision as part of a square canonical
inversion matrix.

The frozen **STRUCTURAL PREDICTION** is:

```text
rank(J_can) = 65,
```

with the weakest singular pair concentrated in the collective-lapse /
homogeneous-scale sector.  Rationale: the known internal null and the known
homogeneous boundary compatibility can pair in the full pre-Legendre map, so
fixing `p_pre` may determine the otherwise free Dirichlet lapse.  This is a
prediction, not an acceptance condition.

If confirmed, the earlier phrase “collective lapse gauge” must be narrowed:
it is a Dirichlet stationary-family degeneracy, not an exact canonical gauge
at fixed initial phase-space data.  Curvature and/or dust would have supplied
a pseudo-constraint that selects it.

## 6. Frozen interpretation of every rank outcome

The calculation must report all 65 singular values and their calibrated
directional errors.  It must not choose a rank from one relative threshold.

- **65 resolved nonzero singular values:**
  `CANONICAL_LEGENDRE_REGULAR`.  The map is locally unique on this reduced
  carrier; an exact canonical lapse null is refuted at the published point.
- **Exactly one error-consistent zero, whose right vector is the analytic
  collective internal lapse with zero `q_new` component:**
  `ONE_CANONICAL_LAPSE_NULL`.  The next geometry may still be unique modulo
  one internal gauge parameter.
- **Additional error-consistent zeros with nonzero final-boundary
  projection:**
  `ADDITIONAL_CANONICAL_DEGENERACY`.  The corresponding boundary directions
  are not uniquely evolved at linear order.
- **Small but resolved nonzero singular values:**
  `RESOLVED_PSEUDOCONSTRAINTS`.  Report the full scale, overlaps and condition
  number; do not round them to gauge.
- **Any weakest direction whose singular value is not separated from its
  calibrated error:**
  `CANONICAL_RANK_NUMERICALLY_OPEN`.
- **Failure of action-gradient reciprocity, Lorentzian branch or derivative
  calibration:** implementation/control failure; no rank claim.

An additional null at the symmetric background is not automatically gauge.
It is labelled `SYMMETRY_OR_GAUGE_OPEN` until compared with a genuinely
stationary nonsymmetric neighbor and supplied with an explicit geometric
generator.  An arbitrary off-shell perturbation cannot settle that question.

## 7. Degrees-of-freedom claim boundary

This census measures local regularity of the reduced Legendre inversion.  It
does not, by itself, prove a continuum physical degree-of-freedom count.

- If rank 65, all thirty reduced boundary configuration directions (one
  homogeneous scale plus 29 schedule-invariant shapes) have a locally unique
  canonical image for admissible initial data.
- If a null vector has final-boundary projection rank `r>0`, at least `r`
  reduced configuration combinations are not uniquely predicted at linear
  order.
- Internal-only nulls may be gauge or history degeneracy and must not be
  subtracted from boundary degrees of freedom without a generator test.
- Modes outside the order-24 invariant quotient, the full 720-edge boundary,
  transverse-traceless graviton identification and the continuum constraint
  algebra remain **NOT TESTED**.

## 8. Framing attack and next boundary

A full-rank Jacobian would establish only a local implicit map.  It would not
show that a finite step exists, that the step is stable, that the scale grows,
or that perturbations follow gravitational-wave dispersion.  Conversely a
rank deficiency at the symmetric point would not establish gauge.

Only after this census may a separate protocol:

1. solve with `p_target=p_pre(published)` and reproduce the published slab;
2. then solve with `p_target=p_post(published)` as the candidate forward
   continuation;
3. recompute the spectrum at that new **on-shell** nonsymmetric point;
4. decompose perturbations into scale and shape modes;
5. later compare shape propagation with a preregistered tensor operator and
   perform a genuine refinement sequence.

No result in the present census supports inflation, a physical clock, `c`, a
Planck scale or particle masses.
