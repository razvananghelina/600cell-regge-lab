# Preregistration: complete-action Lyapunov--Schmidt reduced scalar

Date: 2026-08-13

Prior-art gate and post-result update: `8be130e`

Frozen stagnation result: `1d66278`

Status: **frozen before evaluating a complete-action row at any corrected
state, taking any complete-action Newton step, or inspecting the reduced
scalar after transverse correction**.

The positive scalar seen at the 48 uncorrected anchors is an inspected target.
Any confirmation has **PATTERN-informed provenance**, even though this
protocol freezes its look-elsewhere count and numerical method.

## 1. Complete hypotheses and carrier

Retain without alteration:

- the De Felice--Fabri time-symmetric 600-cell dust sandwich;
- both five-stage schedule parities;
- the order-24 invariant `35=30+5` internal logarithmic carrier;
- the exact collective tangent and its committed `35 x 34` orthonormal
  complement `Q`;
- the four ordered zero-sum Helmert phase contrasts, both signs, and
  `eta=1e-4`;
- all sixteen final-boundary logarithmic vectors;
- the collective grid `t=-0.10,-0.05,0,+0.05,+0.10`;
- the final `z` of each of the 80 states written by result `8f651cd`.

The old boundary, dust mass, external `tau`, simplicial carrier and branch
choices are fixed.  No amplitude reduction, direction selection, restart or
favorable parity is allowed.

For the complete internal logarithmic equation

```text
E = u*(partial S_total/partial u)/24,
F = Q^T E,
g = w(t)^T E,
```

solve `F=0` at fixed `(t,boundary)` and only then classify `g`.  At a
transverse root, `g` is the derivative of the reduced action because the
implicit `z(t)` contribution is proportional to `F` and vanishes.

## 2. Two disjoint complete-action derivative windows

All action calculations use 100 decimal digits.  The nonlinear solver uses
the frozen logarithmic steps

```text
6e-5, 3e-5, 1.5e-5,
```

while independent final validation uses the disjoint steps

```text
2e-5, 1e-5, 5e-6.
```

For either triple `(h1,h2,h3)` and every internal coordinate, reconstruct
the per-edge logarithmic equation from complete-action central differences:

```text
D(h)  = [S(u_i exp(h))-S(u_i exp(-h))]/(48 h),
R12   = (4 D(h2)-D(h1))/3,
R23   = (4 D(h3)-D(h2))/3,
E6    = (16 R23-R12)/15,
delta = E6-R23.
```

Record the full row and error row.  All 210 geometries in an action-row
evaluation must pass the certified Lorentzian branch gates and maximum
imaginary contamination must be below `1e-80`.

The two windows share no step.  Solver values cannot be reused as validation
values.

## 3. Local Jacobian used only as a preconditioner

At every accepted iterate recompute the binary analytic local Jacobian of
`F` in all 34 `z` coordinates using the already frozen steps

```text
5e-4, 2.5e-4, 1.25e-4
```

and the same sixth-order/Richardson construction as protocol `e53dcaf`.
Let the two symmetrized estimates be `J6` and `J23`, with
`delta_J=J6-J23`.

For the complete-action solver residual `F_s` form

```text
p6  = J6^-1 F_s,
p23 = J23^-1 F_s.
```

The preconditioned step is usable only if:

- all 204 displaced Jacobian geometries pass the branch gates;
- both smallest singular values exceed absolute `1e-12`;
- `norm(p6-p23)/max(norm(p6),1e-30)<=0.1`;
- `norm(delta_J p6)/max(norm(F_s),1e-30)<=0.1`.

The Jacobian is never used to classify the physical equation as zero.  Its
positive definiteness is recorded but is not made a gate.

## 4. Frozen complete-action transverse iteration

Start each state at its stored `z`.  Evaluate the solver-window action row.
Project its error as `epsilon_F=norm(Q^T delta)`.  Propagate that error through
the usable local Jacobian:

```text
epsilon_p = norm(J6^-1 Q^T delta).
```

The current transverse system is `SOLVER_TRANSVERSE_ZERO_CONSISTENT` only if

```text
norm(F_s) <= 10*max(epsilon_F,1e-30),
norm(p6)  <= 10*max(epsilon_p,1e-30),
epsilon_p < 1e-5.
```

Otherwise try `z-alpha*p6` at damping factors

```text
alpha=1,1/2,1/4,...,1/1024
```

in that order.  For every trial recompute the complete-action solver row.  A
finite branch-valid trial is accepted if it is itself solver-zero-consistent
after its local Jacobian is recomputed, or if

```text
norm(F_trial)+10 epsilon_F_trial
    < norm(F_current)-10 epsilon_F_current.
```

For the descent comparison, recomputing the trial Jacobian may be deferred
until after acceptance; it is mandatory before a zero label or next step.

Stop at the first of:

- solver transverse zero consistency;
- no accepted damping;
- six accepted iterations.

Record every action row, error, Jacobian certificate, damping and branch
margin.  No SciPy optimizer, finite-difference action Jacobian, Broyden
update, favorable restart or relaxed threshold is allowed.

## 5. Independent validation and reduced scalar

For every solver-zero state, recompute the complete-action row with the
validation window.  With

```text
F_v       = Q^T E_v,
delta_F_v = Q^T delta_v,
p_v       = J6^-1 F_v,
delta_p_v = J6^-1 delta_F_v,
g_v       = w(t)^T E_v,
delta_g_v = w(t)^T delta_v,
```

require all of:

- `norm(F_v)<=10*max(norm(delta_F_v),1e-30)`;
- `norm(p_v)<=10*max(norm(delta_p_v),1e-30)`;
- `norm(delta_p_v)<1e-5`;
- all 210 validation branches pass;
- maximum imaginary contamination is below `1e-80`;
- the binary analytic equation agrees with `E_v` below norm `1e-8`.

Only then is the state `TRANSVERSE_ACTION_VALIDATED`.

For a validated transverse state assign exactly one scalar label:

- `REDUCED_SCALAR_ZERO_CONSISTENT` if
  `abs(g_v)<=10*max(abs(delta_g_v),1e-30)`;
- `REDUCED_SCALAR_RESOLVED_NONZERO` if
  `abs(g_v)>100*max(abs(delta_g_v),1e-30)`;
- `REDUCED_SCALAR_UNRESOLVED` otherwise.

Record the signed scalar, signed error, ratio and sign.  A scalar zero is a
stationary grid hit only after all transverse validation gates pass.

## 6. Frozen sign brackets

For each signed case, inspect the four adjacent grid intervals only after
both endpoints are transverse validated and their reduced scalars are
resolved nonzero.  If their signs differ, perform bisection.

At a midpoint initialize `z` by the arithmetic mean of the two endpoint
solutions, run exactly the transverse solver and independent validation
above, and use only the validated reduced scalar.  Stop at:

- a validated scalar zero;
- an unresolved midpoint;
- interval width below `1e-10`;
- 30 bisections.

If the final width is reached without a scalar-zero validation, record
`BRACKET_LOCALIZED_NOT_ZERO_VALIDATED`; do not count it as a hit.  Do not
search intervals whose endpoint signs agree, so the result is explicitly a
five-point grid-and-bracket scan, not a continuous no-root proof.

## 7. Per-case and aggregate outcomes

A signed case receives:

- `NONLINEAR_STATIONARY_CONTINUATION_FOUND` if at least one validated grid or
  bisection point has reduced scalar zero;
- `NO_STATIONARY_POINT_IN_FROZEN_GRID_AND_BRACKETS` if all five grid points
  are transverse validated, every scalar is resolved nonzero, every required
  bracket resolves under section 6, and no hit occurs;
- `ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED` otherwise.

Report hit fractions over 16 signed cases, 8 direction/parity pairs and 4
phase contrasts, exactly as in the preceding protocols.

Use one global label:

- `NONLINEAR_STATIONARY_CONTINUATION_FOUND` if any signed case hits;
- `SIGN_DEFINITE_REDUCED_SCALAR_ON_FROZEN_GRID` if all 80 grid points are
  transverse validated, all 80 scalars are resolved nonzero with one common
  sign, and no case hits;
- `NO_HIT_MIXED_REDUCED_SIGNS_ON_FROZEN_SCAN` if all sixteen cases receive
  the no-point outcome but the 80 grid signs are not common;
- `ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED` otherwise.

The exact attempt count is the 80 grid states plus every mechanically forced
sign-bracket midpoint; report it.  Do not count the three numerical derivative
steps as separate attempts.

## 8. Acceptance, kill and evidence boundaries

A validated scalar-zero state is a **DERIVED COMPUTATIONAL LOCAL** stationary
continuation on this restricted carrier.  Robustness is its preregistered hit
fraction; one isolated hit remains look-elsewhere limited.

A common resolved nonzero sign at all 80 transverse roots is a **DERIVED
COMPUTATIONAL obstruction on the frozen grid** with **PATTERN-informed
provenance**.  It supports the known higher-order pseudo-constraint mechanism
but does not prove sign-definiteness between grid points.  Therefore it kills
only the frozen grid-and-bracket scan, not the continuous interval or the full
theory.

Any unresolved transverse state, scalar or forced bracket keeps the route
**OPEN NUMERICALLY**.  Solver failure is never called root nonexistence.

No outcome establishes amplitude scaling, quadratic dependence, the other 25
boundary directions, a second slab, the full 840-edge carrier, physical time,
a vacuum, inertia, mass, a causal speed limit or Planck units.

## Operational interruption addendum (2026-08-14)

The first execution was externally interrupted after completing and
validating the forty even-parity grid states and while evaluating the first
complete-action correction batch for the odd parity.  It wrote no result
JSON, no reduced-scalar values were printed or inspected, and no aggregate
outcome was available.  The incomplete in-memory calculation is not treated
as evidence.

Before rerunning, the verifier was amended to write an atomic operational
checkpoint after a whole parity has completed both its frozen solver and
independent validation.  The checkpoint is bound to the protocol commit,
input-result commits, grids, derivative windows, damping list and iteration
limit; an incompatible checkpoint is a hard error.  Resumption skips only an
already completed parity and performs the same downstream checks.  The
checkpoint is deleted after a successful result write.  This changes no
mathematical input, action evaluation, solver step, gate, target or outcome
hierarchy.
