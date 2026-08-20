# First precision run: control thresholds were not conditioning-aware

Date: 2026-08-20

## Frozen failure

- protocol commit: `624fc96`;
- registered implementation: `2871252`;
- first failure artifact: `76e5779`;
- artifact SHA-256:
  `23199cf8da5ed4b41d3022174e75e3035e85ddb1af8b2b9ba5aadf03132d2c68`.

The targeted run ended with

```text
FULL_SCALE_STRUT_PRECISION_CONTROL_FAILED
6/8 checks passed.
```

No full suite was run.

## What passed

The known-answer condition-`10^7` calibration separated direct SVD from
binary64 Gram formation.  In both carrier parities:

```text
GESDD/GESVD discrepancy                       0,
GESVD/80-decimal-Gram discrepancy             4.99e-9, 5.73e-9,
all high-precision eigenvalues                 positive,
all decisive preregistered precision criteria true.
```

These observations are not accepted as a result because two required
controls failed.

## What failed

### Frozen binary-Gram reproduction

The direct condition number was reproduced to every serialized digit, but
the reconstructed binary-Gram discrepancy differed from the old run by

```text
6.15e-11 absolute (even),
7.14e-11 absolute (odd).
```

The protocol compared these already-unstable `~0.03` quantities at `1e-10`
*relative*, an absolute demand near `3e-12`.  This was much tighter than the
normal-equation error scale the same protocol was designed to expose.

### Row-order reversal

Reversing all 1560 rows changed the direct weakest singular values by

```text
2.106e-9 relative (even),
3.138e-9 relative (odd).
```

The frozen threshold was `1e-12`.  For `kappa~=1.10e7`, however,

```text
epsilon_binary64*kappa ~= 2.45e-9.
```

Thus the observed differences are exactly on the forward-error scale of a
backward-stable binary64 calculation.  Row reversal changes floating
summation and bidiagonalization paths; bitwise singular-value invariance was
never a valid expectation at this condition number.

## Classification

**DERIVED HARNESS FAILURE.** The two failed thresholds omitted the condition
number.  They do not falsify the carrier, but the first run cannot be
reclassified retroactively.

**OPEN.** The high-precision resolution remains unaccepted until a disclosed
correction protocol passes.  The correction may change only these two
control bounds, must preserve every decisive criterion, and must retain the
first failure artifact literally.

