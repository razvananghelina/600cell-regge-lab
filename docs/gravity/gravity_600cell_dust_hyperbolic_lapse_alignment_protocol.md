# Preregistration: extreme-tangent versus lapse-subspace alignment

Date: 2026-08-17

Prior-art commit: `24d2ce6`.

Status: explicitly confirmatory and frozen before reconstructing any response
matrix `Y`, selecting any extreme invariant vector, or evaluating any new
principal angle.

## 1. Frozen inputs

Require exact SHA-256 values:

```text
full-boundary tangent JSON:
4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5

all tangent midpoints/radii archive:
816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b

full-boundary tangent verifier source:
c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571

full vertex-lapse Schur artifact:
4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349
```

The tangent artifact must pass `19/19` with outcome
`FULL_BOUNDARY_TANGENT_BLIND_CENSUS_CERTIFIED`; its numeric archive must
contain exactly 224 named arrays.  The Schur artifact must be regular at full
rank `120/120` in both parities.

The exact tick, geometry and local derivative hashes inherited by the
tangent verifier remain required.  Use the same 100-digit geometry, 80-digit
Flint balls, two schedule parities, seven deterministic `2T` sectors and four
derivative variants.  No continuum or speed target may be parsed.

## 2. Geometry-selected weak coordinates

For every parity, select the five pole-edge orbit positions from edge type,
before loading a tangent spectrum.  In a minimal sector of dimension `d`,
these give `5d` weak coordinates inside the `65d` domain

```text
U = (delta internal, delta new).
```

All other `60d` coordinates are strong.  Require exact agreement with the
five positions stored by the prior Schur audit.

For each derivative variant reconstruct the complete projected Hessian and

```text
J = [[K_XX,K_XN],[-K_OX,-K_ON]],
R = [[-K_XO,0],[K_OO,I]],
Y = J^-1 R.
```

Every Flint determinant must exclude zero.  The midpoint `Y` maps the
`60d` boundary phase input into `U`.

## 3. Two independent lapse candidates

Order the domain as `(strong,weak)` and write

```text
J_ordered = [[A,B],[C,D]].
```

Construct the canonical weak Schur lift

```text
C_weak = reorder_back([ -A^-1 B ],
                       [     I    ]),
```

using Flint balls.

Independently reconstruct the geometric vertex-lapse columns `G`:

- one on the selected pole edge;
- `-rho/(exp(s)L0^2-rho)` on each incident internal diagonal;
- zero elsewhere.

Project them with the same minimal `2T` basis.  Require rank `5d` for both
`C_weak` and `G`.  Reproduce every previously stored
canonical-versus-geometric projector distance within `2e-8`.  The two
candidates are not merged or refitted.

## 4. Extreme invariant subspaces

For each frozen tangent midpoint let `k=5d`.  Sort eigenvalue moduli without
loading either lapse candidate.  Define

```text
E_plus  = invariant subspace of the k largest moduli,
E_minus = invariant subspace of the k smallest moduli.
```

Use a reordered complex Schur decomposition.  Its threshold is the geometric
mean of moduli `k` and `k+1` in the corresponding ordering.  Require exactly
`k` selected eigenvalues and a boundary gap ratio greater than two under all
four derivative variants.  This gap gate was fixed before evaluating the
new eigenvectors.

As an independent binary64 control, also construct the corresponding direct
right-eigenvector subspace.  Record its projector distance from the ordered
Schur subspace; do not use either representation selectively after seeing
which aligns better.

Transport both Schur subspaces into slab-domain variations:

```text
U_plus  = colspace(Y E_plus),
U_minus = colspace(Y E_minus).
```

Require rank `5d` after transport.

## 5. Principal-angle comparisons and calibration

For subspaces with orthonormal bases `Q1,Q2`, use

```text
distance = ||Q1 Q1* - Q2 Q2*||2
         = sin(theta_max),
overlap  = minimum singular value of Q1* Q2.
```

Evaluate the four preregistered comparisons per sector:

```text
U_plus  versus C_weak,
U_minus versus C_weak,
U_plus  versus G,
U_minus versus G.
```

For each named comparison let `d_op,d_ops,d_val,d_vals` be the four scalar
distances and define

```text
epsilon_step = |d_op-d_ops| + |d_val-d_vals| + |d_op-d_val|,
```

plus:

- the maximum Schur-versus-direct-eigenvector projector discrepancy after
  transport;
- first-order Flint radius bounds divided by the minimum singular values of
  the transported and candidate column matrices;
- `10 eps_machine` times the largest reported matrix condition involved.

Call their sum `epsilon_distance`.  Assign mechanically

```text
IDENTIFIED       if d_op <= 10 epsilon_distance,
SEPARATED        if d_op > 100 epsilon_distance,
NUMERICALLY_OPEN otherwise.
```

Also print the angle in degrees and the overlap.  A small but resolved
nonzero angle is `SEPARATED` for the exact claim and may separately be noted
as a **PATTERN**; it is not relabelled by visual preference.

## 6. Look-elsewhere accounting

There are exactly

```text
2 parities x 7 sectors x 2 branches x 2 candidates = 56 comparisons.
```

Report all 56 labels and these four fixed hit fractions, each out of 14:

```text
plus/canonical, minus/canonical, plus/geometric, minus/geometric.
```

Do not select the smaller of plus/minus per sector and call that a hit.  A
global identification requires the same branch in all fourteen parity-sector
cases.

## 7. Mechanical outcome

Assign exactly one:

1. `HYPERBOLIC_LAPSE_ALIGNMENT_CONTROL_FAILED` if provenance, geometry,
   response, rank or reproduction controls fail;
2. `HYPERBOLIC_EXTREME_SUBSPACE_OPEN` if any fixed-count Schur selection or
   gap gate fails;
3. `HYPERBOLIC_GEOMETRIC_LAPSE_IDENTIFIED` if one fixed branch is
   `IDENTIFIED` with both canonical and geometric candidates in all 14 cases;
4. `HYPERBOLIC_CANONICAL_WEAK_IDENTIFIED_ONLY` if one fixed branch is
   `IDENTIFIED` with the canonical lift in all 14 cases but the geometric
   condition above fails;
5. `HYPERBOLIC_LAPSE_ALIGNMENT_REFUTED` if all 56 comparisons are
   `SEPARATED`;
6. `HYPERBOLIC_LAPSE_ALIGNMENT_MIXED_OR_OPEN` otherwise.

The verifier passes when it reconstructs and classifies the object.  Only
outcome 3 licenses the statement that the extreme modes are the frozen
geometric lapse space.  Outcome 4 identifies the algebraic weak Schur sector
but leaves its geometric lapse interpretation open.  Outcomes 5 or 6 forbid
dismissing the 119 strong pairs as lapse artifacts.

## 8. Claim boundary

This is a finite coordinate-subspace test, not a gauge theorem.  Even an
exact hit does not show that curvature observables vanish on the modes, and
even separation does not show that they are gravitons.  The next physical
gate would be their action on derived deficit-angle/curvature observables,
followed by a second dynamically solved slab and refinement.

Only the new targeted verifier is run.  The full suite is excluded.
