# Protocol: bounded positive-root search for the refined H4 slab

Date: 2026-08-20

Prior-art gate commit: `29162db`.

This protocol is frozen before evaluating any refined action away from the
inherited fill.

## 1. Frozen inputs

Hash and use exactly:

```text
reproducible/verify_gravity_600cell_refined_h4_stationary_fill.py
  89aab727792e20a81e7577e0425f8fa4b1e84e2a7ae66caa9e79a4aebf3581e7
reproducible/gravity_600cell_refined_h4_stationary_fill.json
  283be37bc7530a3cc4fce9e279272359f107f09fb7b1b0eaff141059bfb4e018
reproducible/verify_gravity_600cell_refined_h4_internal_jacobian.py
  6f74f0a73d15b1e50e61e0afe56d74b162d4a98ac87368979f4bd52fe86b6b4e
reproducible/gravity_600cell_refined_h4_internal_jacobian.json
  b900021c21df67c1de1ae18929be302b0d47d2f267c4a919388711a0a0bf5eaa
docs/gravity/gravity_600cell_refined_h4_internal_jacobian_result.md
  8a6603c810a5615956c8cc2ec8a9b8c3a6e015a0f6495f2f8a58811585418ee5
```

Load only the frozen action definitions from the stationary-fill verifier's
AST and parse the committed combinatorics. Require the upstream outcomes
`REFINED_H4_INDUCED_FILL_OFF_SHELL` and
`REFINED_H4_INTERNAL_JACOBIAN_FULL_RANK_MULTIPLE_CLASSES`, `12` matrix
classes, rank ten and inertia `(8,0,2)` for all 24 schedules.

## 2. Equations and coordinates

For each schedule let `G_sigma(y)` be the ten total-orbit internal log
derivatives of the complete Regge-plus-dust action. Use

```text
y = (log(d_01/d_01^0),...,log(d_23/d_23^0),
     log(rho_0/rho_0^0),...,log(rho_3/rho_3^0)).
```

All physical coordinates therefore remain positive. The six cross entries in
the squared-edge matrix are `+d_rs`; the four vertical entries are
`-rho_r`. The twelve boundary coordinates, mass and `P1` weights are fixed.

For numerical scaling only, use the committed schedule matrix `H_sigma(0)`
and minimize

```text
R_sigma(y) = H_sigma(0)^(-1) Re G_sigma(y).
```

This invertible transformation does not change the zeros. Scientific
acceptance is always checked on `G`, never on the optimizer flag or on `R`
alone.

## 3. Time-reversal reduction control

Before solving, compare every schedule with its reverse at the three frozen
anchors

```text
A0 = 0,
A1 = (0.01,-0.01,0.02,-0.02,0.03,-0.03,-1,-1,-1,-1),
A2 = (-0.03,0.02,-0.01,0.01,-0.02,0.03,-0.5,-1,-1.5,-2).
```

At 80 decimal digits require maximum differences in action and all ten
equations below `1e-60`, with identical branch diagnostics. Only then solve
one representative of each committed class. Otherwise stop with
`TIME_REVERSAL_REDUCTION_FAILED`; do not silently expand or shrink the search.

## 4. Frozen bounded search and look-elsewhere count

The main box is

```text
-0.35 <= y_cross <= +0.35,
-8 <= y_rho <= +2.
```

For each of the 12 representatives use exactly six seeds:

```text
S0 = 0,
S1 = 0.5 p_sigma,
S2 = p_sigma,
S3 = (0,0,0,0,0,0,-4,-4,-4,-4),
S4 = (+.05,-.05,+.05,-.05,+.05,-.05,-1.5,-2,-2.5,-3),
S5 = -S4 with its four rho entries replaced by (-3,-2.5,-2,-1.5),
```

where `p_sigma` is the committed unapplied Newton proposal. Clip a seed only
by `1e-6` inside a bound and record every clipping. This is exactly
`12*6=72` main attempts.

Use SciPy `least_squares` with

```text
method='trf', jac='3-point', diff_step=1e-5,
xtol=ftol=gtol=1e-12, max_nfev=1200, x_scale=1.
```

Evaluate the action at 50 decimal digits. Nonfinite or branch-invalid trial
points receive a fixed finite penalty and are counted. Report every endpoint,
raw and preconditioned residual norm, optimizer status, active bounds,
minimum angle-argument modulus, maximum angle-identity residual and maximum
imaginary contamination. Deduplicate endpoints only for reporting, never for
the denominator of the hit fraction.

## 5. Zero-lapse boundary ladder

Independently of the main outcome, run one additional `S1` attempt for every
class in each box with lower lapse bound

```text
L = -4,-8,-12,-16
```

and the same cross and upper lapse bounds. These are exactly `12*4=48`
attempts. Do not warm-start one bound from another.

Label a class `ZERO_LAPSE_BOUNDARY_PATTERN` only when no finite root is
accepted, every ladder endpoint is active on at least one lower lapse bound,
and its best preconditioned residual decreases strictly as `L` moves from
`-4` to `-16`. This is a **PATTERN**, never a no-root theorem. Report the four
pairwise log slopes; do not fit an exponent.

The complete look-elsewhere denominator is therefore `120` attempts.

## 6. High-precision candidate refinement and validation

An endpoint enters refinement only if it is at least `1e-5` inside every box
face, is branch-valid, and has `norm(R)<1e-7`.

At 100 decimals apply at most 40 damped Newton iterations. Construct the
Jacobian of `G` by centered differences at log step `1e-18`; try damping
`1,1/2,...,1/1024`; accept only strict decrease in `norm(H(0)^-1 G)`. Stop at
`norm(H(0)^-1 G)<1e-35`. No iteration may leave the main box or the branch.

Validate every refined candidate at 140 decimals by:

1. direct analytic `G`;
2. an independent Richardson derivative of the complete action in each of
   the ten log coordinates at steps `1e-10` and `5e-11`;
3. a fresh centered Jacobian at step `1e-15` and its numerical rank relative
   to a precision/step envelope;
4. substitution into the time-reversed schedule.

A `FINITE_POSITIVE_ROOT` requires

```text
norm(H(0)^-1 Re G) < 1e-30,
max_abs(G_analytic-G_action_Richardson)/max(1,max_abs(G)) < 1e-25,
max imaginary action/gradient/curvature < 1e-50,
minimum angle-argument modulus > 1e-8,
maximum angle identity residual < 1e-50,
distance from every main-box face > 1e-5,
reverse-schedule equation difference < 1e-40.
```

Report the root Jacobian rank but do not require full rank for existence.
Candidates within `1e-8` of each other in infinity norm are one distinct root.

## 7. Solver classification controls

Run the identical bounded least-squares wrapper on:

```text
F_positive(x)=diag(1,...,10)(x-x_star),
x_star=(0.1,-0.1,0.05,-0.05,0.02,-0.02,0.2,-0.2,0.1,-0.1),
F_negative(x)=(exp(x_0)+1,...,exp(x_9)+1).
```

The first must recover the known interior root to norm `<1e-10`; the second
must not be classified as a root. Corrupting one schedule matrix to rank nine
must also fail the provenance/rank control before any action solve.

## 8. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_H4_STATIONARY_ROOT_CONTROL_FAILED` for any provenance,
   reconstruction, time-reversal, synthetic-solver or validation-control
   failure.
2. `REFINED_H4_FINITE_ROOTS_ALL_CLASSES` if every class has at least one
   validated finite positive root.
3. `REFINED_H4_FINITE_ROOTS_SOME_CLASSES` if at least one but not all classes
   has a validated root.
4. `REFINED_H4_NO_FINITE_ROOT_FOUND_ZERO_LAPSE_PATTERN` if none has a root
   and every class passes the boundary-pattern definition.
5. `REFINED_H4_NO_FINITE_ROOT_FOUND_OTHER` otherwise.

Outcomes 4 and 5 are bounded computational negatives, not nonexistence
proofs. Only outcome 2 licenses an effective boundary Hessian, after a
separate adversarial replication. Outcome 3 is a temporal-canonicity failure
for the frozen search, not a global theorem. No outcome derives an absolute
tick.

## 9. Deliverables

Write a registered verifier, deterministic JSON, result note, identical
targeted rerun and static registry audit. Do not run the full suite. Do not
compute a boundary Hessian, mode spectrum, `c`, `G` or Planck scale.
