# Preregistration: physical equivalence of homothetic staircase schedules

Date: 2026-08-17

Prior-art commit: `e4428bf`.

Status: frozen before evaluating any affine-hull rank, signature, physical
simplex volume or schedule-equivalence condition.

## 1. Frozen upstream control

Load and byte-hash

```text
reproducible/gravity_600cell_overlay_metric_compatibility.json
SHA-256 7de877b83b5524a1c86788f207ec205fa1eae799ca66bf62c1ae6b46081bb45e.
```

Require its `9/9` controls, its 148 chambers, the exact gcd `r-1`, and its
outcome `OVERLAY_INHERITS_STATIC_METRIC_ONLY`.  These data are not targets for
the new test.  They freeze what must be reconciled: inequality under the
identity parametrization.

## 2. Exact physical vertices

Use symbols `R_minus,R_plus,rho>0`,

```text
Delta=R_plus-R_minus,
T^2=rho+Delta^2,
c=phi/2,
phi=(1+sqrt(5))/2.
```

Represent the spatial vertex directions only through their exact Gram matrix

```text
U_ii=1,
U_ij=c for i!=j.
```

The eight physical vertices in ambient signature `(4,1)` are

```text
p_i^-=(R_minus u_i,0),
p_i^+=(R_plus u_i,T).
```

No numerical scale ratio, lapse, dust mass or continuum trajectory is used.

## 3. Affine hull and Lorentzian signature

First derive

```text
n=(sum_i u_i)/(1+3c),
n.u_i=1,
n.n=4/(1+3c),
```

and require all eight vertices to satisfy exactly

```text
n.x_space-(Delta/T)x_time=R_minus.
```

Independently construct the `8 x 6` homogeneous-coordinate vertex matrix in
the abstract basis `(u_0,...,u_3,time,constant)` and require affine rank four,
equivalently homogeneous rank five, with a one-dimensional affine-hyperplane
nullspace reproducing the formula above.

Require the normal-square identity

```text
N.N=4/(1+3c)-Delta^2/(rho+Delta^2)>0
```

for every `rho>0`.  The inequality must be certified by the exact decomposition

```text
N.N = [4/(1+3c)-1] + rho/(rho+Delta^2),
4/(1+3c)-1 = 3(1-c)/(1+3c)>0.
```

As an independent signature control, introduce affine coordinates

```text
z=x_time/T,
y_i=coefficient of u_i,
sum_i y_i=R_minus+Delta z,
```

eliminate `y_3`, and compute the exact `4 x 4` induced Gram matrix `G_Q` in
coordinates `(y_0,y_1,y_2,z)`.  Its leading spatial block must be
`(1-c)(I_3+J_3)`, and its Schur complement must be

```text
-rho-(3/4)(1-c)Delta^2 < 0.
```

This proves signature `(3,1)` without a floating eigenvalue threshold.

## 4. Exact frustum polytope

Define

```text
Q={0<=z<=1,
   y_i>=0 for i=0,...,3,
   sum_i y_i=R_minus+Delta z}.
```

Enumerate all intersections of four of its six facet hyperplanes and retain
those satisfying all inequalities symbolically for arbitrary positive scales.
Require exactly the eight vertices

```text
(R_minus e_i,0), (R_plus e_i,1), i=0,...,3.
```

Require facet vertex counts `(4,4,6,6,6,6)`, the face-lattice incidence of
`Delta^3 x I`, and four-dimensionality.

The exact coordinate four-volume is frozen analytically as

```text
Vol(Q)=(R_minus^3+R_minus^2 R_plus
        +R_minus R_plus^2+R_plus^3)/24.
```

It comes from integrating the tetrahedral section volume
`(R_minus+Delta z)^3/6`; it is not a fitted target.

## 5. All 24 staircase triangulations

For every order `o=(v_0,...,v_3)` and split `k`, form

```text
conv(p^-_(v_0),...,p^-_(v_k),
     p^+_(v_k),...,p^+_(v_3)).
```

In the reduced `Q` coordinates compute its exact homogeneous determinant and
compare it with the determinant of the corresponding standard parameter
simplex.  Require the orientation ratio

```text
R_minus^k R_plus^(3-k)>0
```

for all 96 labelled simplices.  Hence every physical simplex is nondegenerate
and coherently oriented for every positive pair of scales.

For each schedule require:

- three internal tetrahedral facets, each with incidence two;
- fourteen boundary tetrahedra, distributed as one on each time tetrahedron
  and three on each of the four lateral triangular frusta;
- the sum of its four unsigned simplex volumes equal exactly `Vol(Q)`.

These checks are necessary but are supplemented by an explicit global
non-overlap proof.  Define the projective map from the standard prism with
coordinates `(lambda,tau)` by

```text
d=(1-tau)/R_minus+tau/R_plus,
y_i=lambda_i/d,
z=(tau/R_plus)/d.
```

Require its exact inverse

```text
lambda_i=y_i/(R_minus+Delta z),
tau=z R_plus/[R_minus(1-z)+R_plus z].
```

For positive scales both denominators are positive.  Verify both compositions,
all eight vertex images and the six facet images exactly.  A projective
homeomorphism maps every standard staircase triangulation to a nonoverlapping
triangulation of `Q`; this closes the folding gap left by volume alone.

## 6. Isometry rather than coordinate equality

For every labelled simplex also reconstruct its affine map `A_(o,k)` from
standard parameter coordinates to the reduced physical `Q` coordinates.
Require invertibility.  For every ordered pair of the 32 distinct affine map
types, define the local change of coordinates

```text
H=B^-1 A.
```

Check the exact identity

```text
H^T (B^T G_Q B) H = A^T G_Q A.
```

This is the local tensor law underlying
`h_(o,o')=X_o'^{-1} o X_o`.  Separately require that every `X_o` has the same
restriction

```text
bottom: lambda -> R_minus sum_i lambda_i u_i,
top:    lambda -> R_plus  sum_i lambda_i u_i.
```

Therefore the induced global time-boundary map is pointwise fixed.  No claim
about arbitrary anisotropic data is permitted.

## 7. Edge-length reconciliation

From the common ambient geometry independently rederive:

```text
bottom spatial edge^2 = L_minus^2,
top spatial edge^2    = L_plus^2,
same-vertex strut^2   = -rho,
cross diagonal^2      = L_minus L_plus-rho,
R_minus=phi L_minus,
R_plus =phi L_plus.
```

Require exact equality with the homothetic length formulas already used by the
repository.  This prevents proving equivalence for a different physical
geometry.

## 8. Mechanical outcomes

- All affine-hull, Lorentzian, polytope, triangulation, projective, isometry and
  edge-length controls pass:
  `HOMOTHETIC_SCHEDULES_ONE_LORENTZIAN_FRUSTUM`.
- Coplanarity/signature pass but at least one schedule folds, overlaps, has the
  wrong boundary or wrong volume:
  `HOMOTHETIC_FRUSTUM_TRIANGULATION_FAILED`.
- Any affine-rank, signature, projective-map, isometry, source or exact-
  arithmetic control fails:
  `HOMOTHETIC_FRUSTUM_CONTROL_FAILED`.

The positive outcome **REFUTES the physical interpretation**, not the
arithmetic, of `OVERLAY_INHERITS_STATIC_METRIC_ONLY`: `r-1` remains the exact
obstruction to equality under the identity parametrization, while the 24
metrics are isometric descriptions of one flat frustum for every positive
scale pair.

## 9. Exclusions

Do not evaluate a Regge or dust action, claim anisotropic equivalence, fit an
interior map, infer a clock or run the full verifier suite.  Run only the new
registered verifier.
