# Universal-overlay metric compatibility: result

Date: 2026-08-17

## 1. Provenance

- prior-art gate: `5443238`;
- frozen protocol: `5735a15`;
- pre-evaluation zero-polynomial clarification: `0d2b8b4`;
- registered implementation before first evaluation: `de9a043`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_overlay_metric_compatibility.py`;
- result artifact:
  `reproducible/gravity_600cell_overlay_metric_compatibility.json`;
- artifact SHA-256:
  `7de877b83b5524a1c86788f207ec205fa1eae799ca66bf62c1ae6b46081bb45e`.

Only the targeted verifier was run.  It returns **9/9**.  The full suite was
not run, following the active instruction.

## 2. Mechanical verdict

```text
OVERLAY_INHERITS_STATIC_METRIC_ONLY
```

**DERIVED NEGATIVE.**  The universal 148-chamber overlay is a common
combinatorial refinement of all 24 staircases, but it is not a common metric
refinement for any genuine homothetic evolution.  The complete exact positive
compatibility set is

```text
R_plus/R_minus = 1.
```

Thus the inherited metrics agree only when the two spatial slices have the
same scale.  Neither expansion nor contraction survives schedule independence.

This is not an embedding-only rejection.  The verifier compares the complete
intrinsic Lorentzian pullback Gram matrix on each open overlay chamber.

## 3. Exact construction and controls

For each of the 96 labelled `(order,split)` four-simplices, the verifier solved
the `5 x 5` affine interpolation problem from its five vertices.  It
independently rederived the same map from staircase barycentric weights.  All
96 vertex systems are nonsingular, with determinants

```text
-1: 48,
+1: 48.
```

All `24 x 3=72` internal staircase facets pass exact PL continuity.  The
pullback metric is

```text
g = J^T diag-block(U,-1) J,
U_ii=1,
U_ij=phi/2  (i != j).
```

Every metric is exactly symmetric.  The outer-time/lapse symbol `T` cancels
from every inter-schedule metric difference.  Therefore changing the lapse
cannot restore metric compatibility within this inherited construction.

The static substitution `R_plus=R_minus` sends all 96 affine maps to one
global map and all metrics to one global metric.  The full chamber census is
invariant under the certified `S4 x C2` action.

## 4. Complete look-elsewhere census

There are 32 distinct affine map types among the 96 labelled simplices, and
32 distinct intrinsic metric types.  In this family no two distinct affine
maps accidentally induce the same Gram tensor.

Across the 148 chambers, the number of distinct metrics among the 24
schedules is

```text
 4 distinct metrics :  2 chambers
 6 distinct metrics :  8 chambers
 8 distinct metrics : 28 chambers
10 distinct metrics : 64 chambers
12 distinct metrics : 46 chambers
```

Consequently

```text
identically compatible chambers: 0 / 148.
```

Among the 552 ordered pairs of distinct schedules, their chamber-agreement
counts are

```text
 0 chambers : 336 ordered pairs
19 chambers : 144 ordered pairs
38 chambers :  24 ordered pairs
74 chambers :  48 ordered pairs.
```

No pair agrees everywhere.  This is not a rare bad chamber inside an otherwise
common metric.

## 5. The exhaustive result has a one-line analytic witness

The lexicographically first chamber has sign word

```text
--------------
```

and lies in split `k=3` for both orders

```text
(0,1,2,3),
(0,1,3,2).
```

Write `c=u_i.u_j=phi/2` for distinct tetrahedron vertices and
`Delta R=R_plus-R_minus`.  On these two simplices the spatial affine maps are

```text
X_3 = R_minus sum_i lambda_i u_i + Delta R t u_3,
X_2 = R_minus sum_i lambda_i u_i + Delta R t u_2.
```

Their mixed metric component differs by

```text
(g_3-g_2)_(lambda_0,t)
  = R_minus (R_minus-R_plus) (1-c)
  = R_minus (R_minus-R_plus) (3-sqrt(5))/4.
```

Because `R_minus>0` and `(3-sqrt(5))/4>0`, equality is possible exactly when
`R_plus=R_minus`.  The common time derivative is the same in both maps, so
`T` cancels visibly.

This single exact witness already proves the dynamic no-go.  Independently,
the verifier collected all 19 distinct nonzero compatibility polynomials from
all `148 x 24` assignments.  Both a direct gcd and an intersection of their
irreducible factors give

```text
gcd = r-1,  r=R_plus/R_minus,
```

and direct substitution verifies its only positive root `r=1` in every Gram
matrix.  Hence there is no exceptional nonstatic ratio hidden by a generic
symbolic comparison.

## 6. What this means for gravity

**DERIVED.**  A common subdivision preserves the 24 old PL metrics as 24
different fields on the finer carrier.  It does not identify them.  Therefore
the earlier universal overlay solves the arbitrary **combinatorial** choice of
staircase but does not solve the arbitrary **metric/dynamical** choice.

The proposed route

```text
old staircase metrics -> common overlay -> unique Regge action
```

is closed for every expanding or contracting homothetic slab.  Evaluating a
Regge action on one inherited metric would silently reintroduce the schedule
choice that the overlay was built to remove.

This also explains why the static published sandwich was deceptively benign:
at equal scale the eight boundary coordinates extend to one affine map of the
whole prism, so every triangulation is merely subdividing the same flat local
geometry.  As soon as the scale changes, the top-minus-bottom displacement is
vertex-dependent and the affine interpolants separate.

## 7. Scope: what is not killed

The result does **not** rule out:

- Regge gravity on a deliberately chosen staircase;
- a metric selected directly on the product cell before triangulation;
- a curved-prism or isoparametric element;
- dynamical coarse graining and a perfect action;
- a new variational principle that makes the 32 local metric types dynamical
  and then eliminates them.

All of these add structure beyond common refinement.  None may be called
derived from the overlay alone.

## 8. Post-result primary-source audit

The learned distinction is standard in broad form.  A Regge triangulation
carries a piecewise-affine metric determined simplex by simplex; see
Khatsymovsky, [*Affine connection form of Regge
calculus*](https://arxiv.org/abs/1509.04974), and Mikovic,
[*Piecewise Flat Metrics and Quantum
Gravity*](https://arxiv.org/abs/2001.11439).  Dittrich and Steinhaus explicitly
find triangulation dependence in four-dimensional classical Regge calculus:
<https://arxiv.org/abs/1110.6866>.  Bahr and Dittrich's improved/perfect action
program confirms that discretization independence is a dynamical
coarse-graining problem, not a consequence of a common mesh:
<https://arxiv.org/abs/0907.4323>.

No located source prints the present exact 600-cell tetrahedral-prism witness,
the 148-chamber distribution or the gcd `r-1`.  External novelty remains
**OPEN**; a targeted search cannot prove novelty.

## 9. Next falsification route

The smallest explicit extra principle is the schedule-free tensor-product
interpolant on the prism,

```text
X(lambda,t)
  = ((1-t)R_minus+t R_plus) sum_i lambda_i u_i,
X_time(lambda,t)=T t.
```

It is the unique map that is affine in `lambda` for fixed `t` and affine in
`t` for fixed `lambda`, matches both boundaries and is `S4`-equivariant.  That
conditional uniqueness is promising, but the multi-affinity hypothesis is
new **STRUCTURAL** input and the resulting prism metric is not piecewise flat.

Before using it, the next mission must therefore:

1. audit tensor-product/isoparametric prism elements and curved-simplex Regge
   prior art;
2. prove or refute the stated conditional uniqueness;
3. state explicitly whether affine outer time, proper time or dust time is the
   interpolation coordinate;
4. only then map the certified order-complex vertices and test exact
   Lorentzian nondegeneracy of its fine four-simplices.

If that metric requires another arbitrary time interpolation or produces
degenerate/non-Lorentzian fine simplices, this repair route also closes.
