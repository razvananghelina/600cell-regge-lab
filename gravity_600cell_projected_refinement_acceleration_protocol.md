# Protocol: direct weak-lapse acceleration on projected 600-cell refinements

Date: 2026-08-17

Prior-art commit: `81db1ec`

Status: frozen before evaluating any level-1 or level-2 acceleration
coefficient.

## 1. Question and provenance split

Use the complete direct cellular Regge action fixed in the prior-art note and
ask whether its homogeneous closed-dust acceleration improves under two
projected spatial refinements.

The calculation has two commits.

### Stage A: coefficient production

A registered verifier constructs every carrier and writes all raw weak-lapse
coefficients to

```text
reproducible/gravity_600cell_projected_refinement_acceleration_blind.json
```

It may use the exact level-zero coefficient as an implementation calibration,
but it must not compute distances of refined coefficients from `-1/2`, rank
refinement variants by that distance, or issue a continuum-convergence
verdict.  The artifact and its hash are committed first.

### Stage B: disclosed comparison

Only after the Stage-A artifact commit may a separate registered verifier
read its frozen hash and compare all refined coefficients with `-1/2`.  This
is artifact-level preregistration, not cognitive blindness: the continuum
target and the known coarse excess are already known.

## 2. Frozen carriers and the diagonal ambiguity

The base coordinates and adjacency come from `commons.build_600cell()`.  Its
stored decimal coordinates are first normalized vertex by vertex to unit
norm (the maximum source correction is of order `1e-11`).  The coarse
tetrahedra are the 600 sorted four-cliques.  Shared edge midpoints are then
normalized radially to the unit `S3`.  Thus every strut can consistently have
the same `rho`; retaining unequal source norms would make that hypothesis
false away from a static slice.

The old red-refinement code exposed a genuine regulator ambiguity.  In an
exact regular parent tetrahedron, all three opposite-midpoint diagonals of the
central octahedron have equal length.  The stored decimal coordinates perturb
those equal lengths by about `1e-11`, so raw floating-point minimization is not
a geometric selector.  A single tower is therefore insufficient.

Stage A must construct exactly these four disclosed variants:

1. `legacy_float_shortest`: at both levels use the old raw shortest-length,
   then lexicographic tuple ordering (the algorithm is inherited, while the
   mandatory initial unit normalization means this is not promised to be a
   bitwise copy of the older spectral artifact);
2. `first_tie_rank_0`: in every level-zero parent choose the lexicographically
   first of the three central diagonals, then use raw shortest at level two;
3. `first_tie_rank_1`: choose the lexicographically middle diagonal at the
   first refinement, then raw shortest;
4. `first_tie_rank_2`: choose the lexicographically last diagonal at the first
   refinement, then raw shortest.

The candidate ordering is by the sorted global endpoint pair.  No coefficient
or continuum distance participates in a choice.  All variants must have

```text
(V,T) = (120,600), (840,4800), (6480,38400)
```

and every triangular face must have incidence two.

This four-variant census does not exhaust the `3^600` local resolutions.  It
is a hostile finite regulator control.  Agreement cannot prove universality;
disagreement can already refute a claimed unique finite answer.

## 3. Volume radius, curvature and dust

For a unit-sphere carrier let

```text
Vbar = sum_tetrahedra EuclideanVolume(t),
s0   = (2*pi^2/Vbar)^(1/3).
```

The scaled slice `s0*p_v` has volume radius `R0=1`.  Compute its spatial Regge
curvature independently from Euclidean tetrahedral dihedral angles:

```text
C = sum_edges (s0*l_e) * (2*pi-sum_incident_t theta_t,e).
```

Freeze the conserved total dust mass separately on every carrier as

```text
M=C/(8*pi).
```

The level-zero controls are

```text
s0*edge_length = zeta,
zeta^3 = pi^2*sqrt(2)/50,
C = 720*zeta*epsilon3,
epsilon3 = 2*pi-5*acos(1/3).
```

No local refined dust distribution is claimed.

## 4. Direct irregular cellular action

For each ordered tetrahedron `(0,1,2,3)` with unit positions `p_i`, use the
intrinsic basis

```text
E_i = s_minus*(p_i-p_3),                       i=0,1,2,
E_3 = ((s_plus-s_minus)*p_3, T),
T^2 = rho+(s_plus-s_minus)^2.
```

Its Lorentzian Gram matrix is therefore

```text
g_ij = E_i dot E_j,
g_i3 = s_minus*(s_plus-s_minus)*(p_i-p_3) dot p_3,
g_33 = -rho.
```

Use facet conormals

```text
n_0=(1,0,0,0), n_1=(0,1,0,0), n_2=(0,0,1,0),
n_3=(-1,-1,-1,s_plus/s_minus-1),
n_bottom=(0,0,0,1), n_top=-n_bottom.
```

For conormals `n,m`, the internal-angle cosine is

```text
-<n,m>_(g^-1) / sqrt(<n,n><m,m>).
```

Use the already certified Lorentzian boundary branch

```text
sqrt(-x-i0)=-i*sqrt(x).
```

For a spatial edge of unit chord length `d`, its timelike trapezoid has
signed area

```text
i*(s_minus+s_plus)*d/2
 * sqrt(rho+(s_plus-s_minus)^2*d^2/4).
```

Sum, without angular averaging,

```text
S_grav = -i [
  sum_edges A_edge*(2*pi-sum_tetra theta_lateral)
  + sum_faces A_face_minus*(pi-sum_two_tetra theta_bottom)
  + sum_faces A_face_plus *(pi-sum_two_tetra theta_top)
],

S_total=S_grav-8*pi*M*sqrt(rho).
```

The implementation may use the algebraically identical incidence-expanded
sum `2*pi*sum_unique - sum_local` for vectorization.

## 5. Action controls fixed before refined coefficients

All of the following are mandatory.

1. On level zero, at the three points

   ```text
   (L_minus,L_plus,rho)=(1,1,1/4),(1,4/5,1/10),(1,6/5,1/10),
   s=phi*L,
   ```

   the irregular implementation agrees with the exact closed regular action
   below `5e-10` relative and has imaginary contamination below `5e-9`.
2. At every carrier, static lower/top boundary angles equal `pi/2` within
   `5e-9`, and the static cellular action agrees with `tau*C` below `5e-9`
   relative.
3. With the selected mass, the static total action and logarithmic lapse
   derivative vanish within the errors implied by the derivative audit.
4. Reordering local tetrahedron vertices by a frozen deterministic
   permutation leaves a nonstatic action invariant within `5e-9` relative.
5. Every real evaluation used for a residual has Lorentzian inertia `(3,1)` in
   every tetrahedral frustum.

A failed action or branch control makes the refined coefficient **OPEN**; it
cannot be rescued by a favorable continuum comparison.

## 6. Frozen weak-lapse estimator

Put `eta=tau/R0=tau` and keep `rho=eta^2` for the leading coefficient.  For a
candidate `a`, set

```text
s_minus=s0,
s_plus =s0*exp(a*eta^2).
```

Let

```text
P_minus=(1/2)*partial S/partial log(s_minus),
P_plus =(1/2)*partial S/partial log(s_plus),
F       =partial S/partial log(rho).
```

The seam equation is

```text
P_plus(s0,s0,eta^2)+P_minus(s0,s_plus,eta^2)=0.
```

and the lapse equation is `F=0`.  The omitted `O(eta^2)` lapse correction can
change finite residuals but not their common leading root; disagreement in
that limit is a failure.

The leading seam equation is affine in `a`, rather than an assumed nonlinear
fit.  Indeed, `s_plus-s_minus=O(eta^2)`, timelike areas are `O(eta)`, lateral
angle changes are `O(eta^2)`, and boundary-angle changes are `O(eta)`.  After
the static `O(eta)` dust/curvature cancellation, the one-slab action is
quadratic in the endpoint log-scales through `O(eta^3)`.  Differentiation
with respect to the seam scale therefore gives an affine leading residual.
The different lapse derivative is treated by the preserved correction below.
Stage A must verify the seam affineness numerically at held-out disclosed
values; it may not merely assume it.

Seam logarithmic derivatives are computed by the fixed analytic-branch
complex difference

```text
D_h f = [f(z+i*h)-f(z-i*h)]/(2*i*h),
D      = [4*D_(h/2)-D_h]/3,
h      = 2e-5.
```

At each of the fixed values

```text
a=0,-1,-2,-3,
```

evaluate the primary scaled seam residual `g_eta(a)=G/eta` at

```text
eta = 0.04, 0.02, 0.01, 0.005.
```

For its values `y_0,...,y_3`, define, without regression,

```text
r_i=(4*y_(i+1)-y_i)/3,
q_i=(16*r_(i+1)-r_i)/15.
```

Use `q_1` as the estimate of the limiting affine residual `g_0(a)` and
`abs(q_1-q_0)` as its truncation audit.  With

```text
alpha=g_0(0),
beta =g_0(0)-g_0(-1),
```

the primary coefficient is the exact affine root `a_G=-alpha/beta`.
Require `beta!=0` and verify both held-out affine identities

```text
g_0(-2)-2*g_0(-1)+g_0(0)=0
g_0(-3)-2*g_0(-2)+g_0(-1)=0
```

to relative `2e-6`.  The coarse `q_0` residuals give an independent
truncation coefficient through the same two-point formula.

### Preserved first-run correction: the lapse equation is quadratic

The registered first run in `f474463` correctly passed the affine seam
calibration but falsified the protocol's claim that the lapse equation is
also affine.  Its artifact was preserved in commit `0a57607` before this
correction.

At fixed endpoint difference the lapse derivative contains
`(s_plus-s_minus)^2/rho`.  Since the numerator is `O(a^2*eta^4)` and
`rho=eta^2`, the leading scaled lapse residual is

```text
f_0(a)=kappa_2*a^2+kappa_1*a.
```

The exact static identity forces the root `a=0`; the other root is the
dynamic branch.  Repeat the residual extrapolation using
`eta=0.04,0.02,0.01`, then determine, without regression,

```text
kappa_2=[f_0(-2)-2*f_0(-1)]/2,
kappa_1=kappa_2-f_0(-1),
a_F=-kappa_1/kappa_2.
```

Require `kappa_2!=0`, require the normalized static residual `f_0(0)` below
`2e-6`, and hold out `a=-3` by checking

```text
f_0(-3)=9*kappa_2-3*kappa_1
```

to relative `2e-6`.

The same first run exposed conditioning, not a branch ambiguity: after the
exact `O(eta)` curvature/dust cancellation, the complex action derivative
loses useful digits in `F/eta^3`.  Calibrated only on the known level-zero
answer, the frozen lapse derivative is instead the symmetric real-log rule

```text
D_h^R f=[f(z+h)-f(z-h)]/(2*h),
D^R=[4*D_(h/2)^R-D_h^R]/3,
h_F=2e-3,
```

with a complete repeat at `h_F=1e-3`.  The seam continues to use the original
complex analytic-branch derivative unchanged.  No root is chosen by
proximity to `-1/2`, and there is no regression or adjustable fit window.

The Stage-A coefficient gates are:

- the coefficient from fine versus coarse seam extrapolation differs by less
  than `2e-6`;
- a complete repeat with derivative base step `h=1e-5` differs by less than
  `2e-6`;
- the dynamic lapse root from real derivative steps `2e-3` and `1e-3`
  differs by less than `5e-5`;
- the extrapolated lapse and seam coefficients differ by less than `5e-5`;
- level zero agrees with the already exact radius coefficient
  `-0.5394897340206755...` within `5e-6`.

These tolerances were calibrated only on the known level-zero answer before
any refined coefficient was evaluated.

## 7. Frozen Stage-A outcome

### `PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_DERIVED`

Use **DERIVED NUMERICAL** only if every carrier, action, branch, affine-root,
extrapolation and calibration gate passes.  Print all variants and levels
without a continuum ranking.

### `PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_OPEN`

Use **OPEN** if any definition, branch, numerical stability, uniqueness or
calibration gate fails.  Preserve all diagnostics; do not compare a failed
coefficient with the continuum target.

## 8. Frozen Stage-B comparison and outcomes

For every coefficient `a_(variant,level)` define only in Stage B

```text
ratio=-2*a,
error=abs(ratio-1).
```

Level zero is common.  The finite two-level improvement gate requires, for
all four variants independently,

```text
error(level1) < error(level0),
error(level2) < error(level1).
```

Let `spread_j=max_variant(a_j)-min_variant(a_j)`.  Regulator ambiguity is
called subdominant at level `j` only when

```text
spread_j < min_variant abs(a_j+1/2).
```

This compares choice spread with the smallest remaining discretization
error; it is not a fitted percentage threshold.

The possible disclosed verdicts are:

- `PROJECTED_REGGE_ACCELERATION_ROBUST_TWO_LEVEL_IMPROVEMENT`: every tower
  improves at both steps and the ambiguity is subdominant at levels 1 and 2;
- `PROJECTED_REGGE_ACCELERATION_CHOICE_SENSITIVE_IMPROVEMENT`: every tower
  improves, but the registered regulator spread is not subdominant;
- `PROJECTED_REGGE_ACCELERATION_NOT_UNIVERSAL`: at least one registered tower
  fails strict two-step improvement;
- `PROJECTED_REGGE_ACCELERATION_COMPARISON_OPEN`: the frozen Stage-A artifact
  is invalid or a comparison cannot be made.

Even the strongest verdict is a **DERIVED finite numerical trend**, not an
infinite convergence theorem.  Calling it continuum convergence without
further levels and regulator families remains **PATTERN/OPEN**.

## 9. Scope boundary

This mission tests only the homogeneous Friedmann scale mode with a global
lapse and total dust mass.  It does not derive the tensor Hessian, remove
gauge modes, identify gravitons, determine an effective limiting speed,
select a fundamental tick, or produce Planck or particle scales.  Those
questions remain **OPEN** regardless of the outcome.
