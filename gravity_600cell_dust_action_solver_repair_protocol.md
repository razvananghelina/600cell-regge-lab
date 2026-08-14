# Preregistration: arbitrary-precision complete-action solver repair

Date: 2026-08-14

Prior-art gate: `6b7f9e4`

Frozen numerical-boundary result: `64a13f6`

Status: **frozen before evaluating any new complete-action derivative at a
deformed target state after result `64a13f6`**.  All information from that
result, including the positive raw scalar pattern, is treated as inspected.

## 1. Question, object, carrier and hypotheses

The question is unchanged.  On each of the same 80 frozen
`(parity,direction,sign,t,boundary)` states, solve

```text
F(t,z,b) = Q^T E_action(t,z,b) = 0,        z in R^34,
```

where

```text
E_action,i = (1/24) d S_total / d log(u_i),   i=1,...,35.
```

Only after an independently validated transverse root may the collective
component

```text
g(t,z,b) = w(t)^T E_action(t,z,b)
```

be classified.  `S_total` is the complete Lorentzian Regge-plus-dust action,
not the binary analytic equation.

Retain without alteration:

- the De Felice--Fabri time-symmetric 600-cell dust sandwich;
- both five-stage schedule parities;
- the order-24 invariant `35=30+5` logarithmic carrier;
- the committed `35 x 34` quotient basis `Q` and collective tangent `w`;
- the four ordered zero-sum Helmert contrasts, both signs and `eta=1e-4`;
- all sixteen boundary logarithmic vectors;
- `t=(-0.10,-0.05,0,+0.05,+0.10)`;
- the final `z` and boundary vector written for every state by `64a13f6`.

Every committed binary64 decimal is interpreted as that exact decimal in the
arbitrary-precision calculation.  This repairs numerical differentiation; it
does not silently replace the frozen carrier by a new exact carrier.

The binary analytic local Jacobian `J_b` remains only a proposal generator and
invertible norm preconditioner.  It is not the physical residual or a claimed
Jacobian of the complete action.

## 2. Why the old error proxy is retired

For central differences with expansion

```text
D(h) = D(0) + a h^2 + b h^4 + ...,
```

the old reported value `E6` was sixth order, but its proxy `E6-R23` was fourth
order.  Therefore the observed solver/validation proxy ratio near
`80.98 ~= 3^4` was built into the estimator.  It was not a measured
uncertainty in `E6` and cannot remain a zero or descent gate.

No old scalar or root label is promoted.  Result `64a13f6` remains
`ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED`.

## 3. Frozen high-precision action derivatives

Use 100 decimal digits throughout construction of the variables, complete
action evaluation, central differences, projection by `Q`, and contraction
with `w`.  Do not convert a 35-row to binary64 before those contractions.

For a logarithmic step `h`, define directly

```text
D_i(h) = [S(...,u_i exp(h),...)-S(...,u_i exp(-h),...)]/(48 h).
```

The operational pair is

```text
h_op_primary = 1e-20,
h_op_shadow  = 1e-15.
```

The disjoint final-validation pair is

```text
h_val_primary = 3e-20,
h_val_shadow  = 3e-15.
```

For either pair, use the primary row as `E` and the signed stability proxy

```text
delta_E = D(h_primary)-D(h_shadow).
```

There is no Richardson extrapolation.  At 100 digits the formal primary
central-difference truncation is order `1e-40`, while the shadow difference
exposes an order-`1e-30` term.  These are asymptotic expectations, not rigorous
bounds; every result must therefore retain the label **DERIVED COMPUTATIONAL**.

The two pairs have comparable formal accuracy and share no action point.  The
validation is not made artificially 81 times stricter than the solver.

For every row, independently audit the unperturbed point and all 70 coordinate
perturbations at logarithmic step `1e-6` with the already certified branch
implementation.  Require all negative-count, minimum-Gram and
minimum-angle-argument gates.  Require maximum imaginary contamination of the
high-precision derivative rows below `1e-70`.

## 4. Calibration before any target row

Before evaluating a deformed target state, run both derivative pairs on the
published symmetric control for each parity.

Require all of:

1. the imported published-dust verifier still passes `14/14`;
2. every control branch audit passes;
3. both new primary log-derivative rows agree with the previously certified
   60-decimal action-only control row `u_i*gradient_i/24` in maximum absolute
   norm below `1e-10`; this tolerance covers the old control's frozen
   `3e-6` second-order difference scale and is not a target-root tolerance;
4. operational and validation primary rows agree componentwise within ten
   times the sum of their measured shadow differences, with a `1e-60`
   arithmetic floor;
5. all new imaginary contaminations are below `1e-70`.

If any calibration gate fails, stop before target evaluation with
`DERIVATIVE_CALIBRATION_FAILED`.  Do not change a step or tolerance in the
same experiment.

## 5. Frozen transverse iteration

Start each target at the `final_z` stored by `64a13f6`, represented at 100
decimal digits.  At every accepted iterate:

1. evaluate the operational complete-action row;
2. recompute the frozen binary analytic `J6,J23,delta_J` at the binary64 image
   of that iterate with steps `(5e-4,2.5e-4,1.25e-4)`;
3. form

```text
p       = J6^-1 F,
delta_p = J6^-1 Q^T delta_E,
p23     = J23^-1 F.
```

The preconditioner is usable only under the old frozen conditions:

- all 204 Jacobian branch geometries pass;
- both minimum singular values exceed `1e-12`;
- `norm(p-p23)/max(norm(p),1e-60) <= 0.1`;
- `norm(delta_J p)/max(norm(F),1e-60) <= 0.1`.

Stop as numerically unresolved, without a damping sweep, if the operational
action row fails a branch/imaginary gate, the preconditioner is unusable, or
`norm(delta_p) >= 1e-5`.  These are accuracy failures, not evidence that no
root exists.

An operational row is zero-consistent only if

```text
norm(F) <= 10 max(norm(Q^T delta_E),1e-60),
norm(p) <= 10 max(norm(delta_p),1e-60),
norm(delta_p) < 1e-5,
```

and all row and preconditioner gates pass.

Otherwise choose one active merit before proposing a step:

```text
if norm(p) > 10 max(norm(delta_p),1e-60):
    M = norm(p),       epsilon_M = norm(delta_p)
else:
    M = norm(F),       epsilon_M = norm(delta_F).
```

In the second branch the iterate is not zero-consistent only because the raw
`F` gate failed, so `M-10 epsilon_M` is positive.  This explicit switch also
prevents an algebraic dead zone when only one of the two zero gates passes.

Propose `z-alpha*p` for

```text
alpha = 1,1/2,1/4,...,1/1024
```

in that order.  During a damping sweep keep the current `J6` and the selected
merit fixed.  For the natural merit use

```text
p_trial       = J6^-1 F_trial,
delta_p_trial = J6^-1 delta_F_trial.
```

and for the raw merit use `M_trial=norm(F_trial)` and
`epsilon_M_trial=norm(delta_F_trial)`.  Accept the first branch-valid trial
satisfying

```text
M_trial + 10 epsilon_M_trial < M - 10 epsilon_M.
```

Because an iterate satisfying the right-side uncertainty band is classified
before the sweep, this inequality has no algebraic dead zone at a legitimate
operational zero.  Recompute `J_b` after every accepted step.  Stop at:

- operational zero consistency;
- an action-row or propagated-accuracy failure;
- unusable preconditioner;
- no accepted damping;
- twelve accepted iterations.

No SciPy optimizer, complete-action finite-difference Jacobian, Broyden
update, restart, altered initial state, threshold relaxation or favorable
state selection is allowed.

## 6. Independent final validation

For every operational-zero state, evaluate the disjoint validation pair at
the same high-precision `z`.  Using the final `J6`, require:

```text
norm(F_val) <= 10 max(norm(delta_F_val),1e-60),
norm(J6^-1 F_val) <= 10 max(norm(J6^-1 delta_F_val),1e-60),
norm(J6^-1 delta_F_val) < 1e-5.
```

Also require:

- every validation branch gate and imaginary gate;
- operational and validation primary transverse rows agree within ten times
  the sum of their measured proxy norms, with floor `1e-60`;
- their complete 35-component primary rows obey the analogous gate;
- the old binary analytic equation is recorded only as a diagnostic, never as
  an absolute pass/fail gate at root scale.

The last change is forced by result `64a13f6`: its old `1e-8` gate was much
larger than the residual and compared a full 35-row difference with a
transverse norm.  It certified neither accuracy nor failure.

Only a state passing all gates is `TRANSVERSE_ACTION_VALIDATED`.

For such a state classify the validation scalar exactly as:

- `REDUCED_SCALAR_ZERO_CONSISTENT` if
  `abs(g)<=10 max(abs(delta_g),1e-60)`;
- `REDUCED_SCALAR_RESOLVED_NONZERO` if
  `abs(g)>100 max(abs(delta_g),1e-60)`;
- `REDUCED_SCALAR_UNRESOLVED` otherwise.

Record its sign.  No scalar is classified at an unvalidated state.

## 7. Grid, brackets and outcome hierarchy

Retain the exact sixteen signed cases and 80 grid attempts.  A sign-changing
adjacent interval is eligible for bisection only when both endpoints are
transverse validated and their scalars are resolved nonzero.  If forced, use
the previous frozen midpoint rules, the solver above, at most 30 bisections,
and width floor `1e-10`.  Report every midpoint as an additional attempt.

Per-case and global outcomes are unchanged from protocol `17f9560`:

- `NONLINEAR_STATIONARY_CONTINUATION_FOUND` only for a validated scalar zero;
- `NO_STATIONARY_POINT_IN_FROZEN_GRID_AND_BRACKETS` only for a fully resolved
  signed case without a hit;
- `SIGN_DEFINITE_REDUCED_SCALAR_ON_FROZEN_GRID` only if all 80 roots validate,
  all scalars resolve nonzero with one common sign, and no case hits;
- `NO_HIT_MIXED_REDUCED_SIGNS_ON_FROZEN_SCAN` only under the previously frozen
  complete mixed-sign conditions;
- `ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED` otherwise.

Report hit fractions over 16 signed cases, 8 direction/parity pairs and 4
phase contrasts.  A checkpoint is written atomically after each complete
five-state signed case and is bound to all commits, steps, thresholds, grids,
damping factors and iteration limits above.  It changes no mathematics.

## 8. Acceptance, kill and evidence labels

A validated scalar-zero is a **DERIVED COMPUTATIONAL LOCAL** stationary
continuation on the frozen order-24 carrier.  Because a positive sign was
already inspected, a common nonzero sign remains **PATTERN-informed** even
after this preregistration.

A fully validated common nonzero sign is a **DERIVED COMPUTATIONAL obstruction
on the frozen grid**, not a continuous no-root theorem.  Any unresolved root,
scalar, calibration, descent or forced bracket keeps the route **OPEN
NUMERICALLY**.  Numerical solver failure is never root nonexistence.

No outcome establishes a second slab, full 840-edge stability, a continuum
limit, physical time, vacuum selection, inertia, mass, a causal speed limit,
Planck units or particle masses.
