# Stationary-base audit: first frozen result and coordinate correction

Date: 2026-08-13

Frozen protocol: `61b7f35` (initial commit `fa001af`)

Implementation: `75a3de5`

Registered verifier:
`reproducible/verify_gravity_600cell_dust_stationary_base.py`

Machine-readable result:
`reproducible/gravity_600cell_dust_stationary_base.json`

Targeted run: **12/12 implementation checks passed**.  The full suite was not
run.

## 1. Frozen outcome

The overall frozen verdict is

```text
OPEN NUMERICALLY OR PARITY-SPLIT.
```

The mechanical parity labels were:

```text
even: MULTIPLE_WEAK_SCALE_STATIONARY_BASES_IN_SCAN,
odd : UNIQUE_WEAK_SCALE_STATIONARY_BASE_IN_SCAN.
```

The odd word `UNIQUE` must not be interpreted physically: twelve odd grid
solves were unresolved, so the frozen scan could not exclude further roots.
The protocol assigned the label from the number of passing candidates before
requiring complete localization.  The overall open verdict prevents a false
combined claim, but that precedence is a protocol defect.

## 2. What the run actually observed

For the even schedule, all thirteen frozen points in `[-0.1,0.1]` already
passed the transverse gates at `z=0`; all thirteen scalar values were below
`1e-12`, and all thirteen 100-decimal action audits passed the frozen
weak-scale label.

For the odd schedule, the fixed-preconditioner binary64 iteration resolved
only `t=-0.02,z=0`.  That point passed the action-only audit.  Every trial
remained Lorentzian; the other failures were line-search/precision failures,
not branch failures.

This supplies a **PATTERN** of a collective stationary family but does not
establish it, because the action-derivative empirical errors were about
`2.3e-10` to `3.1e-10`, larger than the extrapolated derivative norms of order
`1e-12`.

## 3. Error found by the independent audit

The mission was motivated by a claimed odd base correction of `8.675e-3`.
That number was wrong.  It projected raw equations

```text
(partial S/partial x)/24
```

onto a Hessian computed from logarithmic equations

```text
x*(partial S/partial x)/24.
```

Inserting the missing `x` factors changes the predicted corrections to

```text
even: 2.521e-10,
odd : 9.019e-7.
```

Both are below the frozen `1e-5` weak-scale tolerance.  The alleged large odd
displacement is **RETRACTED / REFUTED BY COORDINATE CONSISTENCY**.  This is the
main result of the round: the hostile action-only check exposed an error in
the hostile framing itself.

## 4. Remaining numerical limitation

The action-only sixth-order estimator used steps

```text
2e-4, 1e-4, 5e-5.
```

At those steps the empirical Richardson difference is still about `2.6e-10`
near `t=0`.  Therefore a derivative of order `1e-12` is consistent with zero
but not resolved at its reported magnitude.  The old and new gates show that
no large correction is required; they do not yet prove stationarity across
the whole interval on a scale below the soft modes.

## 5. Correct next test

Preregister an action-only scan of all thirteen collective points for both
parities with steps

```text
2e-5, 1e-5, 5e-6.
```

At 100 decimal digits there is no roundoff reason to retain the larger steps;
the fourth-order Richardson difference should fall by approximately `1e4`.
Do not run a binary64 transverse solver.  At every point compute the complete
35-component logarithmic gradient, its projection through the frozen
quotient inverse, parity agreement and branch margins.

Only if that smaller-step audit resolves all points may the collective path
be called a stationary family at the weak scale.  Nonlinear boundary
continuation remains downstream.
