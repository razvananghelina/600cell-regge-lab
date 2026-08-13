# Preregistration: weak-scale action-only collective-path audit

Date: 2026-08-13

Prior-art and coordinate correction: `ad5f0ad`

First stationary-base result: `ad5f0ad`

Precision quotient result: `29a779f`

Status: **frozen before evaluating any complete action at the three smaller
steps**.

## 1. Question and hypotheses

On each derived even/odd five-stage schedule, keep the old and final regular
600-cell boundaries, unrounded source `l0`, reconstructed dust mass, and the
35 orbit-reduced internal carrier fixed.  At the thirteen preregistered
collective points

```text
t = -0.10, -0.075, -0.05, -0.03, -0.02, -0.01,
     0,
     0.01, 0.02, 0.03, 0.05, 0.075, 0.10,
rho(t)=tau^2 exp(t),  q(t)=l0^2-rho(t),
```

is the complete 35-component logarithmic action gradient consistent with
zero at a numerical precision fine enough to resolve displacement through
the four `4.605e-8` soft modes?

No transverse root solver, mass fit, schedule selection, boundary
perturbation or comparison target is allowed in this correction.

## 2. Known and open boundary

Lapse freedom and higher-order pseudo-constraints are **KNOWN STRUCTURE** in
Regge calculus.  Exact action flatness along the displayed collective path is
already **DERIVED COMPUTATIONAL**.  The first action-gradient audit used
steps whose empirical error was about `2.6e-10`, so weak-scale stationarity of
the path remains **OPEN**.

The coordinate-consistent base correction estimates `2.521e-10` and
`9.019e-7` are frozen controls, not targets to reproduce.

## 3. Frozen action-only derivative

At every `t`, parity, internal logarithmic coordinate `i`, and step

```text
h1=2e-5, h2=1e-5, h3=5e-6,
```

evaluate the independently certified complete total action at 100 decimal
digits on

```text
x_i -> x_i exp(+h),  x_i -> x_i exp(-h).
```

Define the per-edge logarithmic equation

```text
D_h[i] = [S(+h)-S(-h)]/(48 h),
R12=(4 D_h2-D_h1)/3,
R23=(4 D_h3-D_h2)/3,
E6 =(16 R23-R12)/15,
epsilon=norm(E6-R23).
```

All 5460 action points must remain Lorentzian with minimum Gram modulus above
`1e-8` and minimum angle-argument modulus above `1e-6`.  Maximum action or
derivative imaginary contamination must be below `1e-80`.

## 4. Frozen weak-mode conditioning

Use the already committed exact base Householder complement `Q` and each
schedule's precision-corrected `34 x 34` quotient matrix `H_Q`; do not
rediagonalize or modify them.  Let

```text
e       = E6,
de      = E6-R23,
delta   = norm(H_Q^-1 Q^T e),
epsilon_delta = norm(H_Q^-1 Q^T de),
w(t)    = normalized collective path tangent,
g       = w(t)^T e,
epsilon_g = abs(w(t)^T de),
floor_e = max(epsilon,1e-30),
floor_delta=max(epsilon_delta,1e-30),
floor_g=max(epsilon_g,1e-30).
```

The frozen adequacy bounds for a subsequent internal perturbation scale
`1e-4` are

```text
epsilon_delta < 1e-5,
delta         < 1e-5,
norm(e)       < 1e-10,
abs(g)        < 1e-11.
```

## 5. Per-point labels

Assign exactly one:

### STATIONARY_WITHIN_ACTION_ERROR_AND_WEAK_SCALE

All branch/imaginary gates and all adequacy bounds pass, and

```text
norm(e) <= 10 floor_e,
delta   <= 10 floor_delta,
abs(g)  <= 10 floor_g.
```

### RESOLVED_SMALL_NONZERO_GRADIENT

All adequacy bounds pass, but at least one of the three quantities exceeds
`100` times its corresponding floor.  The printed point remains an adequate
linear base but is not called stationary.

### RESOLVED_NONSTATIONARY_ON_WEAK_SCALE

`delta > 1e-5`, `delta > 100 floor_delta`, and the numerical/branch gates
otherwise pass.

### WEAK_SCALE_PATH_NUMERICALLY_UNRESOLVED

Every other case.  No threshold may be loosened after inspection.

## 6. Parity and cross-step controls

At each `t`, compare even and odd `E6` rows.  Report their norm difference
and the combined empirical envelope

```text
epsilon_parity = epsilon_even+epsilon_odd.
```

Label them `PARITY_AGREES_WITHIN_ACTION_ERROR` if the difference is at most
`10 epsilon_parity`, `PARITY_RESOLVED_DIFFERENT` if it exceeds
`100 epsilon_parity`, and unresolved otherwise.  Equality is not presumed.

Also require the median ratio of the old first-run error to the new error at
the shared thirteen even points to exceed `1e3`; otherwise the intended
precision correction is unresolved.  The expected factor `1e4` is reported,
not a PASS target.

## 7. Outcomes and claim boundary

Per parity assign:

- `ALL_13_PATH_POINTS_STATIONARY_WITHIN_ERROR`;
- `ALL_13_PATH_POINTS_WEAK_SCALE_ADEQUATE_SOME_NONZERO`;
- `PATH_HAS_RESOLVED_NONSTATIONARY_POINT`;
- `PATH_AUDIT_NUMERICALLY_UNRESOLVED`.

If both parities receive the first label and all thirteen parity comparisons
agree, the result is **DERIVED COMPUTATIONAL ON THE FROZEN GRID** and a
**PATTERN** of a continuous stationary family.  Thirteen samples do not prove
the entire interval or an analytic gauge identity.

Even the strongest result does not show that collective lapse remains gauge
after boundary deformation.  It does not recertify the boundary mixed block,
establish nonlinear continuation, cover the 840-edge carrier, or select a
clock or Planck scale.
