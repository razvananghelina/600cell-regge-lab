# Canonical-lapse first run: unresolved Jacobian reporting failure

Date: 2026-08-16

Prior-art gate: `c7f3e29`.

Frozen protocol: `ded77c5`.

Implementation-before-evaluation commit: `3a90633`.

## Recorded first run

The targeted verifier returned **7/7** with the mechanical outcome

```text
CANONICAL_LAPSE_JACOBIAN_OPEN.
```

Both parities reproduced the fixed-lapse seed and its uniform canonical
mismatch, then refused to start Newton because the calibrated initial
two-variable Jacobian did not meet the frozen rank-two/error-band gate.

The first-run artifact SHA-256 was

```text
46320e3fd6df6b28e5593a8ef4cbb9f236cc2a4ef38599a21cae884bbf5b4485.
```

## Reporting defect

The Jacobian was evaluated and used for the decision, but the implementation
serialized Jacobians only after an accepted Newton step or at a converged
endpoint.  Since the initial matrix failed before any accepted step, the
artifact contains

```text
history=[];
endpoint_jacobian=null;
failure=JACOBIAN_ERROR_BAND_FAILURE.
```

Consequently the artifact does not expose the singular values, matrix
entries or error estimate needed to distinguish an exact/structural rank-one
signal from a merely unresolved numerical calibration.

## Status and frozen correction

The scientific outcome is **PROVISIONAL / NUMERICALLY OPEN** until the
already computed decision matrix is recorded.  This is not permission to
change a step, tolerance, rank gate, seed, target, action or solver.

The sole correction is to store every attempted Jacobian, including the
first unresolved one, in a separate `jacobian_attempts` list and rerun the
same frozen verifier from the beginning.  No outcome rule changes.
