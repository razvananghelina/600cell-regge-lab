# Recorded numerical failure: published dust control in binary64

Date: 2026-08-13

Frozen protocol commit: `cc0902b`

First targeted run: **8/10 implementation checks passed**.

This file records the failed run before any numerical correction.

## What passed

- the upstream 65-variable audit remained `33/33`;
- the unrounded source formulas reproduced the paper's printed sandwich with
  maximum squared-length discrepancy `7.64e-14`;
- both schedule parities remained nondegenerate Lorentzian;
- the full 2400-simplex and 100-orbit actions agreed to relative
  `3.51e-9` (even) and `3.77e-9` (odd);
- their 65 analytic gradients agreed to relative `4.62e-11` and `4.94e-11`.

The analytic sourced per-edge residuals were:

```text
even: max pole     4.945377440890e-12
      max diagonal 1.509686420411e-14

odd:  max pole     4.514883282051e-10
      max diagonal 2.120421620425e-13
```

Thus both parities met the preregistered `1e-7` physical stationarity gate and
were provisionally classified `FULL_REPRODUCTION`.

## What failed

The binary64 centered differences of the **total** action failed:

```text
even: relative derivative error 1.2706004e-3
      imaginary derivative      9.9329083e-6

odd:  relative derivative error 9.99999999e-1
      imaginary derivative      5.3555098e-2.
```

At the published point

```text
rho=tau^2=0.00010404,
delta=rho*(3e-6)=3.1212e-10.
```

The sourced pole derivative is a cancellation between gravitational and dust
orbit derivatives of order `10^3`, while the desired total derivative is
order `10^-8` or smaller.  The action itself also has a binary64 branch
residual of order `10^-7`; dividing its differenced noise by `delta` explains
the observed amplification.  This is a diagnosis, not permission to discard
the failed preregistered gate.

## Status

- **DERIVED NUMERICAL FAILURE:** the frozen binary64 finite-difference control
  does not certify the total derivatives.
- **PROVISIONAL PATTERN:** all analytic full/reduced gradients satisfy the
  published full-sandwich equations.
- **NOT YET CERTIFIED:** `FULL_REPRODUCTION`.
- **NOT A PHYSICS NEGATIVE:** the failure occurs in the independent numerical
  derivative, not in the stationary residual itself.

The correction must keep the same geometry, source, action, step and physical
threshold and replace only the arithmetic used for the independent action
difference.
