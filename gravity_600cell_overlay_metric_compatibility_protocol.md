# Preregistration: metric compatibility of the universal staircase overlay

Date: 2026-08-17

Prior-art commit: `5443238`.

Status: frozen before computing any dynamic affine map, Gram matrix, chamber
census or compatible scale ratio.

## 1. Frozen inputs

Use only the committed artifacts

```text
reproducible/gravity_600cell_universal_staircase_overlay.json
SHA-256 0dd03eed878f599463a44160484c74ddeaa0511fc70c8b2e77bc05a2f36dd3dc

reproducible/gravity_600cell_overlay_face_poset.json
SHA-256 439a3d067d50415f0a47c79091ec746c12dd7975b2246b6143f3f7a70847ce13
```

and require their recorded outcomes and pass counts.  The first supplies all
148 strict sign words and their unique simplex assignment for each of the 24
orders.  The second independently supplies the exact 33 arrangement vertices,
148 chamber closures and full face-poset controls.

No result count, compatible ratio or example from the new calculation may be
inserted into these inputs.

## 2. Frozen boundary geometry

Keep `R_minus`, `R_plus` and `T` as independent exact symbols during the
affine and symmetry computations, with `R_minus,R_plus>0`.  For the four unit
directions of one regular 600-cell tetrahedron use the exact spatial Gram matrix

```text
U_ii = 1,
U_ij = phi/2 for i != j,
phi = (1+sqrt(5))/2.
```

The bottom and top images in the abstract basis `u_0,...,u_3` are

```text
bottom i : (R_minus e_i, 0),
top i    : (R_plus e_i, T).
```

The inherited five-dimensional embedding has

```text
T^2 = rho+(R_plus-R_minus)^2,
```

but `T` is common to all schedules.  It must be retained while deriving each
metric and then shown mechanically to cancel from every inter-schedule metric
difference.  No numerical value of `R_minus`, `R_plus`, `rho` or `T` is fitted
or sampled.

Use reduced domain coordinates

```text
x=(lambda_0,lambda_1,lambda_2,t),
lambda_3=1-lambda_0-lambda_1-lambda_2.
```

## 3. Independent affine reconstruction

For every order `o=(v_0,v_1,v_2,v_3)` and split `k`, form the five parameter
vertices

```text
bottom v_0,...,bottom v_k,
top v_k,...,top v_3.
```

Primary construction: solve the exact `5 x 5` affine interpolation system for
each of the four abstract spatial coordinates and the time coordinate.  Reject
any singular vertex matrix.

Independent control: derive the map from staircase barycentric weights.  If
`p=v_k`, then on the simplex

```text
alpha_i = lambda_i                         for i before p,
beta_i  = lambda_i                         for i after p,
beta_p  = t-sum_(i after p) lambda_i,
alpha_p = sum_(i at/after p) lambda_i-t,

X_space = R_minus sum_(i before p) alpha_i u_i
          +(R_minus alpha_p+R_plus beta_p) u_p
          +R_plus sum_(i after p) beta_i u_i,
X_time  = T t.
```

Require exact equality between the independently solved affine forms and this
closed form for all 96 labelled `(order,k)` pairs.  Also require equality of
neighboring affine maps on each of the `24 x 3` shared staircase facets.

## 4. Exact pullback metric

For each affine map obtain its constant Jacobian `J`.  In the abstract target
basis use

```text
eta = diag-block(U,-1),
g_(o,k)(R_minus,R_plus,T) = J^T eta J.
```

Every entry is an exact polynomial over
`Q(sqrt(5))[R_minus,R_plus,T]`.  Require symmetry,
and require the coefficient of every positive power of `T` to vanish in every
difference `g_(o,k)-g_(o',k')`.  This checks rather than assumes that lapse and
the homothetic outer-time correction cannot repair a schedule mismatch.

For each of the 148 chambers and each order, independently reconstruct the
unique split directly from its sign word.  Compare all 24 exact Gram matrices
at the same labelled parameter coordinates.  Record:

1. the full multiset of the number of distinct affine maps per chamber;
2. the full multiset of the number of distinct metrics per chamber;
3. the number of chambers on which all 24 metrics agree identically in both
   spatial scales;
4. the agreement count for every ordered pair of schedules;
5. one lexicographically first exact mismatch witness, if one exists.

Affine-map inequality is diagnostic only.  The verdict is based exclusively
on metric equality.

## 5. Static and symmetry controls

Substitute `R_plus=R_minus` into every map and metric.  Require all 96 affine
maps to reduce to the same global map

```text
X_space=R_minus sum_i lambda_i u_i,  X_time=T t,
```

and all chamber/order metrics to agree exactly.  A failure here is a control
failure, not a physical result.

Reconstruct the full `S4 x C2` action on sign words.  For spatial permutations
the Gram matrices transform by the induced exact coordinate Jacobian.  For
time reflection use `t -> 1-t`, swap `R_minus` with `R_plus`, and reverse the
target time coordinate.  At minimum, require the number of distinct schedule
metrics to be constant on every chamber orbit.  Record any stronger equality
controls that can be checked without choosing a schedule.

## 6. All compatible scale ratios

Only at this stage normalize `R_minus=1`, write `R_plus=r`, and require `r>0`.
For every chamber choose the lexicographically first
schedule as reference and collect every nonzero entry of every one of the 23
metric differences.  Regard them as univariate polynomials in `r` over
`Q(sqrt(5))` after the verified cancellation of `T`.

Compute their exact monic greatest common divisor.  Independently factor every
polynomial and intersect its real positive root set.  The two methods must
give the same complete set.  Verify every candidate by direct substitution in
all `148 x 24` Gram matrices.  Roots `r<=0` are inadmissible; `r=1` is the
static control, not a dynamic success.

This global polynomial test is the acceptance boundary.  A generic symbolic
mismatch alone is insufficient because it could miss an exceptional positive
dynamic ratio.

## 7. Mechanical outcomes

- All controls pass and at least one exact compatible root `r>0`, `r!=1`
  survives direct substitution on every chamber and schedule:
  `OVERLAY_INHERITS_DYNAMIC_REGGE_METRIC`.
- All controls pass, the complete positive compatible root set is exactly
  `{1}`, and no chamber is identically compatible across all schedules:
  `OVERLAY_INHERITS_STATIC_METRIC_ONLY`.
- All controls pass, the complete positive compatible root set is exactly
  `{1}`, and between one and 147 chambers are identically compatible across
  all schedules:
  `OVERLAY_METRIC_COMPATIBILITY_PARTIAL_ONLY`.
- Any source, affine reconstruction, continuity, static, polynomial-root or
  symmetry control fails:
  `OVERLAY_METRIC_COMPATIBILITY_CONTROL_FAILED`.

The first outcome is **DERIVED METRIC**, not yet a Regge equation.  The second
or third is a **DERIVED NEGATIVE** for inheritance from the old staircases,
not a no-go theorem for a separately selected symmetric/perfect action.

## 8. Exclusions

Do not average schedules, choose an order, assign barycentric physical
coordinates, construct a Regge/dust action, solve an evolution equation, fit a
continuum metric or run the full verifier suite.  Only the newly registered
targeted verifier may be evaluated in this mission.

## 9. Pre-evaluation zero-polynomial clarification

There is one logically possible case in which Section 6's nonzero polynomial
list is empty: all 24 metrics agree identically on all 148 chambers for every
positive pair `(R_minus,R_plus)`.  In that case define the global gcd to be the
zero polynomial, report the compatible ratio set as `ALL_POSITIVE_RATIOS`, and
return `OVERLAY_INHERITS_DYNAMIC_REGGE_METRIC` after all other controls pass.

This convention is fixed before implementation or evaluation.  It prevents an
unexpected strongest-positive result from being misclassified as a polynomial
control failure.
