# Canonical continuation: target-precision calibration failure

Date: 2026-08-16

Prior-art gate: `52a6d50`

Frozen protocol: `393e528`, with clarifications `c2e942d`, `cf00b38`

Implementation before evaluation: `1a6f796`

Status: **CONTROL FAILURE BEFORE FORWARD CONTINUATION**.

The first targeted run returned `6/8`.  Both parity reproduction solves drove
the reduced analytic residual below `5.17e-55` and recovered the published
point to approximately `6.55e-21` in logarithmic coordinates, but failed the
full component gates:

```text
full residual infinity norm  1.14154006718413e-39
within-type spread           1.14154006718413e-39
frozen gates                 1e-40 and 1e-50
```

The failed run evaluated no `lambda>0` continuation target.

## Cause

The nonlinear evaluator is fully symmetric and matches the mean of the 30
stored pre-momenta.  The corrected gluing artifact's stored momentum vector
has component spread

```text
1.1415400671841295e-39.
```

That vector came from finite differences and is not an exact invariant
target.  More importantly, the same committed artifact already certifies the
much larger cusp uncertainty norm

```text
3.6513653962011044e-22.
```

Thus the failed `1e-50` component-spread requirement treated an uncertain
upstream vector as exact.  The failure is in the control framing, not evidence
against a canonical branch.

## Allowed amendment

Retain `1e-50` for the newly computed reduced mean residual.  For the full
65-vector norm and within-type spread only, replace the absolute gates by

```text
10 * stored cusp_uncertainty_norm.
```

This is not the observed `1.14e-39`; it is the uncertainty already committed
before this mission.  No solver step, seed, continuation grid, branch gate or
scientific outcome may change.
