# Protocol: nested tangential vertex-displacement prolongation

Date: 2026-08-22

Status: frozen before implementation and before evaluating any new
prolongation, finite-difference derivative or induced edge variation.

Prior-art commit: `79b612b`

## 1. Frozen inputs

Use `/home/razvan/science/.venv/bin/python` and freeze these exact inputs:

```text
commons/cell600.py
  840d921355e040bd4125dc8f8688b9702d63d9119e6f955f6e40b444c2d7d7a7

reproducible/verify_gravity_600cell_projected_rank_edgewise_carrier.py
  50876c582cf22d86296f3f2b715ff1cf3276a9f1320baa3b37d365ce91f2aa23

reproducible/gravity_600cell_projected_rank_edgewise_carrier.json
  b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84

reproducible/verify_gravity_600cell_refined_canonical_map_feasibility.py
  36fba835048e6e0f0676b749192a9d882406932770a00ba1396929bbc4d04a32

reproducible/gravity_600cell_refined_canonical_map_feasibility.json
  ab6209bc745b4c988b59b8c0416522dd2e4a434f17f4cfd596df817bb48ff02e
```

The new verifier may import `build_600cell` from `commons` but may not import
or execute an older verifier.  It must reconstruct both spatial carriers.

## 2. Frozen object

Reconstruct

```text
K0 = P(sd K_600),
K1 = P(Esd_2(sd K_600)).
```

For every unordered fine key `(i,j)`, with `i=j` for a retained old vertex,
define

```text
y_(i,i) = x_i,
y_(i,j) = (x_i+x_j)/||x_i+x_j||.
```

For a coarse tangent field `u`, define `T u` by

```text
(T u)_(i,i) = u_i,
(T u)_(i,j) = (I-y y^T)(u_i+u_j)/||x_i+x_j||.
```

The domain is the direct sum of the tangent spaces of the unit `S^3` at all
`2640` coarse vertices, of dimension `7920`.  The codomain is the analogous
fine tangential carrier.  No momentum variable or action Hessian is part of
this protocol.

## 3. Mechanically distinct carrier routes

### Route A: edgewise-facet traversal

Enumerate the eight exact `r=2` edgewise facets of every ranked barycentric
chamber.  Merge their fine vertices by the unordered repeated-vertex key
used in the accepted carrier construction.  Preserve every occurrence of a
key in a parent chamber so parent consistency can be checked.

### Route B: old-plus-edge census

Independently enumerate all edges of `K0` directly from its tetrahedra.  Form
the fine key set as

```text
{(i,i): i in V(K0)} union E(K0).
```

This route may not call the edgewise-facet enumerator.  It reconstructs fine
positions directly from singleton and edge keys.  The unordered key sets and
positions from the two routes must agree after key sorting.

The analytic derivative is checked against direct centered differentiation
of the nonlinear normalized-point construction.  This is the mechanically
different derivative route; it may not call the analytic derivative helper.

## 4. Frozen deterministic probes and tolerances

Use two nonzero deterministic coarse tangent fields:

```text
u_i^(a) = (I-x_i x_i^T) z_i^(a),
z_i^(0) = (1, ((i mod 7)-3)/7, ((i mod 11)-5)/11,
           ((i mod 13)-6)/13),
z_i^(1) = (((i mod 5)-2)/5, 1, ((i mod 17)-8)/17,
           ((i mod 19)-9)/19).
```

Normalize each complete field to unit Euclidean norm.  Use centered steps
`epsilon in {2^-18,2^-20}`.  Normalize every perturbed coarse vertex before
forming every perturbed refined point.

Frozen gates:

```text
unit-position residual                    <= 2e-10
coarse/fine tangency residual             <= 2e-10
route-A/route-B coordinate residual       <= 2e-10
parent-occurrence coordinate residual     <= 2e-10
analytic/finite-difference vertex error   <= 2e-7
analytic/finite-difference edge error     <= 3e-7
O(4)-covariance residual                  <= 2e-10
old-vertex restriction residual           = 0 in index arithmetic,
                                           <= 2e-15 numerically
```

The finite-difference error is the maximum absolute component error over all
fine vertices or all fine squared-edge derivatives.  It must decrease when
the step changes from `2^-18` to `2^-20`, unless both values are already below
`2e-10`; this prevents a loose tolerance from being the only control.

Test covariance under both fixed transformations

```text
R1 = signed coordinate cycle (0,1,2,3) -> (1,2,3,0), last sign negative,
R2 = I - 2 w w^T,  w=(1,1,1,1)/2.
```

Both must first satisfy `R^T R=I` within `2e-15`, and one must have
determinant `+1` and the other `-1`.  If the stated signed cycle has the wrong
determinant, change only which single output sign is negative so the two
declared determinant classes are represented; record the actual matrices.

## 5. Exact injectivity certificate

Do not infer injectivity from an SVD.  Route B must certify that every old key
`(i,i)` occurs exactly once.  Define the restriction `R_old` by looking up
these keys.  Check in index arithmetic that

```text
R_old T = identity
```

on all coarse tangent blocks.  This proves rank `7920` without a numerical
threshold.  The output must distinguish this vertex-space rank from the rank
of an induced edge-length carrier, which is not asserted here.

## 6. Fine-edge derivative

Enumerate all `134880` fine spatial edges from the fine tetrahedra.  For each
probe field compare

```text
d ||y_a-y_b||^2 = 2 (y_a-y_b) . ((T u)_a-(T u)_b)
```

with centered differentiation of the complete nonlinear coarse-to-fine
construction.  This checks that the prolongation induces an unambiguous fine
metric variation; it does not claim that the resulting edge carrier is the
full Regge configuration space.

## 7. Parent and schedule independence

For Route A, recompute the position and derivative of every occurrence of a
fine key from its actual parent chamber and require equality after merging.
For temporal schedules, freeze the accepted feasibility artifact and require
at both levels:

```text
schedule_count = 24,
all schedules have the same spatial boundary-edge count,
the spatial key/prolongation construction reads no schedule order.
```

The fact that internal temporal diagonals differ is retained; it is not
silently averaged away.

## 8. Negative controls

Both controls must fail the relevant invariant:

1. replace one midpoint by weights `(2*x_i+x_j)/||2*x_i+x_j||`; Route A and
   Route B must then disagree by more than `1e-4`;
2. delete one old key before forming `R_old`; the exact left-inverse census
   must fail.

The verifier must also assert that it loads no continuum eigenvalue,
dispersion, graviton, `c`, `G`, Planck or particle target and constructs no
action Hessian.

## 9. Frozen checks and outcome hierarchy

The registered verifier must report exactly these fourteen checks:

1. frozen provenance;
2. `K0` topology and unit geometry;
3. Route-A fine census and topology;
4. Route-B fine census;
5. route agreement;
6. nonzero midpoint denominators and parent independence;
7. coarse and fine tangency;
8. exact old-vertex left inverse and rank `7920`;
9. centered vertex-derivative agreement;
10. centered fine-edge derivative agreement;
11. two-class `O(4)` covariance;
12. schedule independence with temporal ambiguity retained;
13. both negative controls discriminate;
14. target firewall and unique outcome.

Outcome hierarchy:

```text
NESTED_VERTEX_PROLONGATION_CONTROL_FAILED
NESTED_VERTEX_PROLONGATION_NOT_CANONICAL
NESTED_VERTEX_PROLONGATION_DERIVATIVE_FAILED
NESTED_TANGENTIAL_VERTEX_CARRIER_DERIVED
```

Only the last outcome advances the feasibility gate.  It must be labelled
**DERIVED COMPUTATIONAL / STRUCTURAL INFRASTRUCTURE**, not a restored
constraint.

## 10. Kill and acceptance boundaries

- **KILL:** failure of key uniqueness, parent independence, tangency,
  covariance, the exact left inverse or the metric derivative closes this
  declared nested-carrier route.  No fitted eigenvector matching may replace
  it.
- **ACCEPT:** all fourteen checks pass and the artifact is byte-identical on
  a second targeted run.  This establishes only one canonical nested
  tangential configuration carrier.
- **STILL OPEN AFTER ACCEPTANCE:** a matched on-shell `K1` finite-height
  background, a dynamically selected cotangent relation, normal/lapse
  displacement transport, constraint restoration, gravitons, `c`, `G`,
  Planck units and particle physics.

Run only this verifier and the static registry/coverage guard required to
register it.  Do not run the full suite.
