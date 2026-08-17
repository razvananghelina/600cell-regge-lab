# Preregistration: universal staircase-overlay chamber census

Date: 2026-08-17

Prior-art commit: `4fb243b`.

Status: frozen before testing feasibility of any of the `2^14` sign patterns.

## 1. Frozen local object

Use the closed prism

```text
P = {lambda_i >= 0, sum_i lambda_i = 1, 0 <= t <= 1}
```

and the 14 affine forms

```text
h_A(lambda,t) = t - sum_(i in A) lambda_i,
empty != A != {0,1,2,3}.
```

The full-dimensional open chambers are the nonempty intersections of the
strict interior of `P` with one choice of `h_A>0` or `h_A<0` for every `A`.
Because each intersection is convex, every feasible strict sign pattern is
exactly one chamber.  This identifies the chamber count without numerical
point sampling or connected-component heuristics.

## 2. Exact exhaustive enumeration

Use exact rational linear arithmetic in Z3.  Starting from

```text
lambda_i > 0, sum_i lambda_i = 1, 0 < t < 1,
```

perform a depth-first binary traversal over all 14 signs.  At each node add
the next strict sign constraint and ask satisfiability.  An unsatisfiable
prefix prunes only its descendants.  Record every feasible depth-14 sign
tuple; this is exhaustive over all `2^14=16384` possibilities.

No floating-point tolerance, fitted coordinate, target chamber count or
random sample is permitted.  Store the sorted feasible sign words in the
result artifact so the census can be independently reproduced.

## 3. Staircase common-refinement test

Exhaust all 24 permutations `v=(v_0,v_1,v_2,v_3)`.  For each feasible chamber,
identify each `k=0,...,3` whose strict interior conditions hold:

```text
sum_(i>k) lambda_(v_i) < t < sum_(i>=k) lambda_(v_i).
```

When the lower tail is empty, its condition is the already imposed `t>0`;
when the upper tail is full, its condition is `t<1`.  Otherwise the lower
condition is the positive sign of its subset form and the upper condition is
the negative sign of its subset form.

Require exactly one `k` for every `(chamber, permutation)` pair.  Record the
four assignment multiplicities for every permutation.  This tests whether
every overlay chamber lies in exactly one open four-simplex of each staircase.

## 4. Symmetry controls

Act by all 24 permutations of the tetrahedron vertices.  If `g` is a
permutation, transport the sign on `A` to the sign on `g(A)`.  Require the
feasible chamber set to be invariant.

For interval reflection `t -> 1-t`, use the exact identity

```text
h_A(lambda,1-t) = -h_(A^c)(lambda,t).
```

Thus transport the sign on `A` to the opposite sign on `A^c`.  Require
invariance.  Enumerate exact chamber orbits under the resulting `S4 x C2`
action and record their sizes.  Every orbit size must divide 48 and the sizes
must sum to the total chamber count.

## 5. Face-restriction and gluing control

For each `j`, restrict symbolically to `lambda_j=0`.  Map each four-vertex
subset `A` to `A\{j}` among the remaining three vertices.  Require that,
after duplicate forms are identified, the restricted internal forms are
exactly

```text
t - sum_(i in B) lambda_i,
empty != B != remaining three vertices,
```

while `B=empty` and `B=all` are only the prism-boundary forms `t=0` and
`t=1`.  Require the same multiset structure for all four faces.  This is the
mechanical compatibility condition needed for local overlays to glue on the
600-cell boundary; it does not enumerate the global face poset.

## 6. Independent analytic controls

Require:

1. exactly 14 distinct nontrivial subset forms and complement pairing without
   fixed points;
2. all 24 vertex permutations and both time parities produce 48 distinct
   transformations of labelled sign words;
3. each staircase order uses exactly the nested tail chain with sizes
   `3,2,1` as its three internal walls;
4. the union of internal walls over all orders is all 14 forms;
5. every recorded sign word is rechecked by a fresh exact solver instance,
   independent of the prefix traversal.

## 7. Mechanical outcomes

- All controls pass, every chamber has one assignment for every order, and
  every face restriction passes:
  `UNIVERSAL_STAIRCASE_OVERLAY_CERTIFIED`.
- Enumeration and base controls pass, but at least one chamber/order has zero
  or multiple assignments:
  `UNIVERSAL_STAIRCASE_OVERLAY_NOT_COMMON_REFINEMENT`.
- Any arithmetic, exhaustive-census, independent recheck, symmetry, orbit or
  face-restriction control fails:
  `UNIVERSAL_STAIRCASE_OVERLAY_CONTROL_FAILED`.

The first outcome is **DERIVED COMBINATORIAL** only.  A large chamber count
must be reported as a practical cost, not hidden behind the positive label.

## 8. Exclusions

Do not construct or evaluate a Lorentzian Regge action, dust action, Hessian,
canonical map, diffusion dimension or continuum fit in this mission.  Do not
run the full verifier suite; run only the newly registered targeted verifier.

