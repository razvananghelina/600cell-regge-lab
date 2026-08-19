# First universal-local-lift execution: reconciliation control failure

Date: 2026-08-19

## Frozen execution state

The protocol was committed as `dd302d8` and the first registered source as
`ac5f80f`.  The first targeted execution stopped during the old-block
reconciliation for the baseline `(lambda,tau)=(2,5)` construction.  No JSON
artifact was written and no local block coefficient or final scientific
verdict was printed.

The traceback was:

```text
RuntimeError: the old radial ansatz does not reproduce the declared data up to metric sign
```

The run had already entered the branch in which the baseline 48-unknown
affine system returned a solution.  Therefore strict blindness to existence
of a baseline block is no longer claimed.  The full block, its support,
residual, other constructions, and final outcome were not disclosed.

## Cause

The verifier incorrectly identified the second 120 canonical-data
coordinates with unit normal-displacement amplitudes `nu`.  The exact global
lift instead uses the raw strut squared-length variations `s`.

The frozen old carrier derives, at each vertex,

```text
edge datum = 8 lambda (sigma_u + sigma_v),
s_v        = 6 (lambda - 1) sigma_v - 2 tau nu_v.
```

Thus the two 240-dimensional data domains are equal for nonzero `tau`, but
their bases differ.  In the new `(sigma,s)` basis,

```text
nu_v = 3 (lambda - 1) sigma_v / tau - s_v / (2 tau).
```

Substituting the raw radial-plus-normal displacement block without this basis
change made the reconciliation control ask the wrong equality.

## Frozen correction before rerun

Do not alter the 48-unknown affine constraints, candidate data map, support
criterion, corruption, global residual, or outcome hierarchy.

Change only reconstruction of the old physical block:

1. build its original `(sigma,nu)` radial-plus-normal displacement matrix;
2. right-compose it with the exact basis change above to `(sigma,s)`;
3. under the reversed overall metric-sign convention only, apply the unique
   overall sign needed for its local Jacobian image to equal the declared
   positive `(sigma,s)` data block;
4. then compare the new and old flex blocks and their Poincare-kernel
   difference exactly as preregistered.

This is a control-coordinate correction, not a scientific-coefficient edit.
Preserve this first failure permanently and freeze the first post-correction
artifact before any further interpretation.

