# Preregistration: precision correction for the canonical-lapse Jacobian

Date: 2026-08-16

Prior-art gate: `c7f3e29`.

Original protocol: `ded77c5`.

Unresolved first result: `d30854c`.

Status: **frozen before evaluating the corrected Jacobian**.

## 1. Reason for the correction

The first run returned `CANONICAL_LAPSE_JACOBIAN_OPEN`.  Its operational
Jacobian had smallest singular value `4.2445844786e-9`, while its calibrated
error was `1.2636262572e-8`.  The determinant printed by that run is therefore
not evidence of rank two.

The same weak scale was independently resolved in the upstream full canonical
Legendre calculation with much smaller, calibrated central-difference steps.
This correction tests whether the open verdict was caused by truncation error.
It does not reinterpret or overwrite that verdict.

## 2. Frozen mathematical problem

Retain without change the exact object, carrier, Regge+dust evaluator,
conserved mass, two schedule parities, sole seed, two unknowns `(s,z)`, reduced
equations `(F0,F1)`, committed post-momentum target, branch gates, Newton rule,
damping order, iteration bound, residual tolerances, complete 65-component
substitution, parity test and outcome hierarchy of `ded77c5`.

In particular, no new seed, coordinate rescaling, mass, target, root bracket,
optimizer, fitting parameter or lapse prior is allowed.

## 3. Sole permitted numerical change

Change only the four coordinate-wise central-difference steps to the values
already used successfully for the weak canonical scale:

```text
operational primary = 1e-20
operational shadow  = 1e-15
validation primary  = 3e-20
validation shadow   = 3e-15
```

Keep 100-decimal arithmetic and the original entrywise calibration and rank
gate:

```text
epsilon = ||d_op||_2 + ||d_val||_2 + ||d_cross||_2 + 1e-60,
s_min > 100*epsilon.
```

The corrected implementation must record all four Jacobians and their error
diagnostics at every attempted Newton state.

## 4. Mechanical outcomes

The original seven-outcome hierarchy remains binding.  The precision run is
interpreted as follows:

- a control or calibration failure leaves the route **OPEN NUMERICALLY**;
- an unresolved smallest singular value leaves the route **OPEN NUMERICALLY**
  and forbids Newton;
- a resolved rank-two Jacobian licenses exactly the already-preregistered
  deterministic Newton solve;
- only `HOMOTHETIC_CANONICAL_LAPSE_SELECTED` after all original gates licenses
  a local two-variable canonical root.

No closeness of the seed residual, nonzero printed determinant, Newton
trajectory or approximate root counts as success.

## 5. Interpretation boundary

A successful run would be **DERIVED COMPUTATIONAL LOCAL** only within the
homothetic two-variable restriction.  It would show that canonical
consistency selects a relative next lapse together with a next scale for the
fixed initial data.  It would not establish an absolute clock, full
65-variable uniqueness, gauge restoration, refinement stability, continuum
gravity or a fundamental time unit.  The lapse remains **STRUCTURAL / candidate
pseudo-constraint** pending those tests.

Only the registered targeted verifier may be run; the full suite is excluded
by the user's instruction.
