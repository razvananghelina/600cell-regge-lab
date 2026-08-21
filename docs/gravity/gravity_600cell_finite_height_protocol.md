# Preregistration: exact finite-height fixed-state root census

Date: 2026-08-21

Prior-art gate commits: `c8edca2`, corrected before calculation by `726e52f`.

Status: frozen before evaluating the elimination function `D(v,q)`, its
factorization, any finite positive root or any numerical root scan.

## 1. Frozen action and state family

Use the accepted primary and adversarial cubic artifacts, without importing
their finite-root conclusions:

```text
gravity_600cell_generic_velocity_cubic.json
SHA-256 1d35b46cd4db20df0af3ed3e6b5de676d69753cf5059e0eb607d1eec949b9103

gravity_600cell_generic_velocity_cubic_adversarial.json
SHA-256 b5167d597a927f8b441a096c31034aa04efa435883284dc2d9bfbd3b9cb3ff0d
```

For every real `v`, including `v=0` and both roots of `K(v^2)`, freeze

```text
L_minus=1,
M=mu(v),
p0=p(v).
```

Solve the exact positive-height equations

```text
F(1,L_plus,h^2;mu(v))=0,
p_pre(1,L_plus,h^2;mu(v))=p(v),
h>0,
L_plus>0.
```

No Taylor endpoint, fitted height, changed state, cosmological constant or
extra degree of freedom is permitted.

## 2. Exact affine reduction before root search

Differentiate the complete action first.  Then introduce

```text
q=(L_plus-1)/h,
L_plus=1+hq,
rho=h^2.
```

Derive the exact residuals and certify

```text
partial_h^2 C=0,
partial_h^2 P=0.
```

Only after that certificate, extract

```text
C=C0(v,q)+h Ch(q),
P=P0(v,q)-p(v)+h Ph(q).
```

Form independently, with the displayed sign convention,

```text
D(v,q)=C0(v,q) Ph(q)-[P0(v,q)-p(v)] Ch(q).
```

Record exact formulas for all four coefficients and `D` before any sampled
evaluation.  Verify directly that `(h,q)=(0,v)` cancels both boundary
residuals for every real `v`.

## 3. Complete algebraic case split

No division by a slope is allowed before the following cases are classified
on the full real `q` line:

1. `Ch!=0`: set `h=-C0/Ch`, require `D=0`, then substitute into both original
   residuals;
2. `Ch=0,C0!=0`: no solution;
3. `Ch=0,C0=0,Ph!=0`: solve the momentum equation and substitute;
4. `Ch=Ph=0`: determine whether constants are compatible, inconsistent or
   leave a continuum of heights.

Every real zero of `Ch`, `Ph` and their gcd/resultant must be enumerated
exactly or certified by rigorous interval bounds.  Squared or
denominator-cleared candidates must be substituted into the unsquared
expressions.

For every surviving root require exactly

```text
h>0,
1+hq>0.
```

Record the Jacobian determinant of `(C,P)` with respect to `(h,q)` at every
root.  A nonzero determinant proves an isolated update for fixed `v`; a zero
determinant requires a separate local multiplicity or continuum analysis.

## 4. Trivial boundary factor and cubic control

Before searching away from `q=v`:

1. verify `D(v,v)=0` exactly;
2. determine the exact multiplicity of the factor `q-v` without assuming it;
3. divide only by the certified multiplicity;
4. show how the first nonzero diagonal expansion coefficient relates to the
   already-certified positive cubic cross-resultant

   ```text
   Delta(v)=129600 epsilon(v)^2/(v^2+4).
   ```

This is a control that the finite equation contains the local no-jet theorem.
It may not be used to infer the sign of the remaining global factor without
proof.

## 5. Global root classification

The root census must cover

```text
(v,q) in R^2
```

including `v=0`, the two `K=0` velocities and both unbounded directions.
A complete result requires one of:

- an exact factorization plus analytic sign/monotonicity proof;
- a theorem reducing all roots to finitely many certified one-dimensional
  equations with complete real-root isolation;
- outward-rounded interval arithmetic on a compactification of every
  remaining domain, together with analytic tail bounds.

Suggested compactification for diagnostics is

```text
v=sinh(s),
q=sinh(t),
```

but a finite `(s,t)` box without tail proofs is not a census.  A numerical
grid, arbitrary-precision root finder or plotted sign map is **PATTERN** only
and forces the global outcome `FINITE_HEIGHT_OPEN` unless upgraded to one of
the complete proofs above.

Before any exploratory scan, print the exact diagonal multiplicity, all
slope exceptional sets and the asymptotic limits of the normalized
elimination function as `q->+-infinity` for fixed `v`.

## 6. Frozen numerical controls

After the exact formulas exist, compare the complete full-action residuals
with their affine reconstructions at 100 decimals for

```text
v in {-3/2,0,2/3,5/2},
q in {-2,-1/3,1/2,3},
h in {1/7,2/5},
```

retaining only triples with `1+hq>0`.  Require normalized differences below
`1e-80`.

Hostile controls:

1. `M -> mu(v)+1/10` must change `C0` by exactly `-4*pi/5`;
2. `p0 -> p(v)+1/10` must change the momentum constant by exactly `-1/10`;
3. reversing the sign in the determinant definition must fail at the frozen
   point `(v,q)=(2/3,-1/3)` unless both terms separately vanish there.

The hostile point is a control chosen before evaluation; it may not be moved
after the first run.

## 7. Outcome hierarchy

### `FINITE_HEIGHT_NO_POSITIVE_ROOTS`

Use **DERIVED NEGATIVE, scoped** only after a complete all-real classification
proves that every non-boundary root fails `h>0`, endpoint positivity or one of
the original equations.  This closes the finite homogeneous tick route for
the fixed action.

### `FINITE_HEIGHT_ISOLATED_UPDATES`

Use **DERIVED EXACT / STRUCTURAL** if a complete classification gives one or
more isolated positive roots.  Report the exact root count as a function of
`v`, every bifurcation value and the full multiset of branches.  Call them
state-dependent discrete updates, not ticks.

### `FINITE_HEIGHT_CONTINUUM`

Use **STRUCTURAL reparametrization** if a positive-dimensional family of
positive roots survives on a nontrivial state domain.

### `FINITE_HEIGHT_EXCEPTIONAL_STRATA`

Use **STRUCTURAL / OPEN physical selection** if roots occur only at isolated
or lower-dimensional state values and the set is completely classified.

### `FINITE_HEIGHT_OPEN`

Use **OPEN** for incomplete real-root isolation, missing tail bounds, failed
controls or disagreement between exact and numerical routes.

No isolated root is a fundamental tick until composition, stability and
carrier/action refinement pass.  Global scale covariance limits every result
to a dimensionless ratio; no outcome here derives seconds, `c`, `G` or Planck
time.  Only the new targeted verifier will be run.
