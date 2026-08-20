# Correction protocol: conditioning-aware precision controls

Date: 2026-08-20

This target-disclosed correction is committed after the `6/8` first failure
and before rerunning the precision verifier.  It changes no matrix, spectrum,
arithmetic precision, decisive comparison or outcome hierarchy.

## Frozen failure

| input | SHA-256 |
|---|---|
| first-failure note | frozen by the commit containing this protocol |
| first-failure artifact | `23199cf8da5ed4b41d3022174e75e3035e85ddb1af8b2b9ba5aadf03132d2c68` |
| original precision protocol | `603aa2bd2c54de143df3598b7d5d03cac07338de51d5b646d068dcef2498d7e2` |

The artifact must retain outcome
`FULL_SCALE_STRUT_PRECISION_CONTROL_FAILED`, `6/8`, zero direct-driver
disagreement, positive high-precision spectra and decisive criteria true in
both parities.

## Correction 1: reconstruction control

Retain the direct condition-number reproduction criterion

```text
relative(kappa_reconstructed,kappa_frozen) < 1e-10.
```

Replace the ill-posed relative comparison of two inaccurate binary-Gram
errors by

```text
abs(error_reconstructed-error_frozen)
    < epsilon_binary64 * kappa_reconstructed^2.
```

This bound is the already-preregistered normal-equation uncertainty scale;
it is not inferred from the `6e-11` observed difference.

## Correction 2: row-order control

Replace the absolute `1e-12` relative threshold by the dimension-aware
backward-stability envelope

```text
row_reversal_relative_error
    < 240 * epsilon_binary64 * kappa_GESVD.
```

The factor 240 is the frozen column dimension, not a fitted safety factor.
The protocol now asks whether row reversal stays within the ordinary
conditioned forward-error scale, not whether two floating algorithms are
bitwise invariant.

## Unchanged gates

All of the following remain byte-for-byte conceptually unchanged:

- known-answer `delta=1e-7` calibration;
- two direct LAPACK drivers;
- exact-binary 80-decimal Gram accumulation and eigensolve;
- `1e-8` direct/high-precision comparison;
- strict positivity and `1e40` precision margin;
- `0.05--20` normal-equation explanation ratio;
- deletion-of-pole corruption;
- outcome hierarchy and interpretation firewall.

If either corrected control fails, the result remains a control failure.  If
the corrected controls pass but a decisive criterion changes, the result is
a precision disagreement.  Only a complete pass may resolve the old
`NUMERICALLY_OPEN` qualifier.

