# Protocol: nested 6+3+1 reduction of the refined H4 equations

Date: 2026-08-20

Prior-art gate commit: `7714933`.

This protocol is frozen before evaluating a new reduced branch.

## 1. Frozen inputs

Hash and use exactly:

```text
reproducible/verify_gravity_600cell_refined_h4_stationary_root.py
  0105508a17fc40f34eb7a15f4c7c36bc89850c653476cb6a2e8086ca4281cceb
reproducible/gravity_600cell_refined_h4_stationary_root.json
  e945dc54a0768b00358aca6bef9e9a105ab3d0080d22dd83dfd140b038adf14d
docs/gravity/gravity_600cell_refined_h4_stationary_root_result.md
  c41c81409e2aa8d16bd7db71e68c3e954e5eaf568e4618ea52153262824b42ff
docs/gravity/gravity_600cell_refined_h4_nested_reduction_prior_art.md
  fe395c93a1ce5209fccc4829e23011be185b3a1ec5c5674c0d68bce11f6a9f0d
```

Load definitions only from the root verifier's AST. Require its exact
`13/13`, `0/120` bounded-negative outcome, 12 schedule classes and all fast /
high-precision cross-evaluator controls. Reuse its frozen action,
combinatorics, branch convention and complex128 evaluator without changing a
formula.

## 2. Exact decomposition

Let `x` be the first six cross-log coordinates and `z` the last four rank-
lapse logs. Define the normalized Helmert matrix

```text
Q = [(1,-1,0,0)/sqrt(2),
     (1,1,-2,0)/sqrt(6),
     (1,1,1,-3)/sqrt(12)].
```

Parameterize

```text
z=t*1_4+Q*u,  u in R^3.
```

For schedule `sigma`, define

```text
F_sigma(x,z)=first six components of G_sigma(x,z),
R_sigma(z)=last four components of G_sigma(x_sigma(z),z),
C_sigma(t,u)=Q^T R_sigma(t*1+Q*u),
g_sigma(t)=1_4^T R_sigma(t*1+Q*u_sigma(t))/2.
```

Here `x_sigma(z)` solves `F=0` and `u_sigma(t)` solves `C=0`. A finite scalar
zero is accepted only by substituting into all ten original equations.

## 3. Frozen linear controls

For every class reconstruct from the committed high-precision Hessian

```text
A=H_xx,
S=H_zz-H_zx A^-1 H_xz,
B=Q^T S Q.
```

Require rank `(6,4,3)` for `(A,S,B)`, the committed cross smallest singular
value `5.93908093...`, `cond(S)<1.1`, inertia `(2,0,2)` for `S`, and equality
of all 12 `S` matrices within their committed error envelopes. A deliberately
rank-deficient `A` must fail this control.

## 4. Inner six-equation solve

At fixed `z`, solve `F(x,z)=0` in

```text
-0.35 <= x_i <= +0.35
```

using the frozen fast evaluator and residual `A^-1 F`. Use SciPy
`least_squares` with

```text
method='trf', jac='3-point', diff_step=1e-5,
xtol=ftol=gtol=1e-12, max_nfev=600, x_scale=1.
```

Try exactly two deterministic seeds: `x=0` and the clipped base linear
predictor `x=-A^-1 H_xz z`. Select the smaller final `norm(A^-1 F)`, counting
both attempts. The inner solve is resolved only if the endpoint is
branch-valid, lies more than `1e-6` inside its box and

```text
norm(A^-1 F)<1e-8.
```

Optimizer success alone does not count.

## 5. Contrast continuation

Use common-lapse grid

```text
t=-7.5,-7,...,+1.5
```

but solve in the fixed continuation order

```text
-2,-1.5,-1,...,+1.5,-2.5,-3,...,-7.5.
```

At every trial require each component of `z=t*1+Q*u` to remain in `[-8,+2]`.
Solve `C(t,u)=0` using residual `B^-1 C` and

```text
method='trf', jac='3-point', diff_step=1e-5,
xtol=ftol=gtol=1e-11, max_nfev=400, x_scale=1,
-4 <= u_j <= +4.
```

At `t=-2` use `u=0`. Thereafter use the nearest resolved continuation point.
If that attempt is unresolved, make exactly one fallback attempt from the
clipped base linear predictor

```text
u_lin(t)=-B^-1 Q^T (R(0)+S(t*1_4)).
```

Do not continue through an unresolved point when constructing brackets. A
contrast point is resolved only when its selected inner solve is resolved,
the full geometry is branch-valid, all four lapse logs and three contrast
coordinates are interior by `1e-6`, and

```text
norm(B^-1 C)<1e-7.
```

## 6. Independent grid validation

At every resolved grid point evaluate the original ten equations at 80
decimals. Require fast/high-precision relative gradient disagreement below
`5e-9`, physical imaginary contamination below `1e-50`, minimum angle
argument above `1e-8`, and high-precision norms

```text
norm(A^-1 F)<1e-6,
norm(B^-1 Q^T G_rho)<1e-6.
```

Record all four `G_rho` components and `g(t)` at high precision. Repeat the
original equations on the time-reversed schedule and require difference below
`1e-40`.

## 7. Scalar candidates

Within each schedule separately, candidates are:

- resolved grid points with `|g|<1e-7`;
- adjacent resolved grid points in sorted `t` with opposite signs and both
  endpoint magnitudes above `1e-7`.

Refine a bracket by at most 30 bisections, solving the nested equations at
each midpoint from the averaged endpoint coordinates. Stop at width `1e-8`
or `|g|<1e-9`. Never bridge an unresolved midpoint.

Every candidate then enters the frozen 100-decimal full ten-equation damped
Newton refinement and 140-decimal action-derivative validation from the root
verifier. Only that existing `FINITE_POSITIVE_ROOT` gate accepts a root.
Deduplicate accepted roots within `1e-8` in infinity norm.

## 8. Controls and outcomes

Run the same nested wrapper on a synthetic system with known cross branch,
known contrast branch and scalar root at `t=-1`; it must recover the root to
`1e-9`. A synthetic scalar `exp(t)+1` must produce no root.

Use the first applicable outcome:

1. `REFINED_H4_NESTED_REDUCTION_CONTROL_FAILED` for provenance, linear,
   evaluator, symmetry, synthetic or validation failure.
2. `REFINED_H4_NESTED_FINITE_ROOTS_ALL_CLASSES` if every class has a validated
   finite root.
3. `REFINED_H4_NESTED_FINITE_ROOTS_SOME_CLASSES` if only some classes do.
4. `REFINED_H4_NESTED_BRANCH_UNRESOLVED` if no root is accepted and any class
   has fewer than all 19 resolved grid points.
5. `REFINED_H4_NESTED_NO_ROOT_GRID_SIGN_CONSISTENT` if all grid points resolve,
   no root is accepted and every class has one nonzero scalar sign.
6. `REFINED_H4_NESTED_NO_ROOT_GRID_OTHER` otherwise.

Outcomes 4--6 are not continuous no-root theorems. Outcome 5 licenses a
separate interval-sign attempt; it is not itself an exclusion certificate.
Any schedule mismatch remains a canonicity problem.

## 9. Deliverables

Write a registered verifier, deterministic JSON, result note, identical
targeted rerun and static registry audit. Do not run the full suite or compute
an effective boundary Hessian, spectrum, tick, `c`, `G` or Planck scale.
