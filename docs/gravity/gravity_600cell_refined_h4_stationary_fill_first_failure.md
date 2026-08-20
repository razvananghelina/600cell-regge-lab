# First failure: refined H4 stationary-fill census

Date: 2026-08-20

Protocol commit: `3d36c54`  
Implementation/registration commit: `b2b3523`

## Scope reached

The first targeted execution passed provenance, exact `K0` topology, exact
rank geometry and conditional `P1` rank-mass controls.  It explicitly built
all 24 schedule slabs and found, for every schedule,

```text
pentachora       = 57,600,
triangles        = 149,280,
triangle types   = 28,
mixed H4 types   = 0.
```

It stopped at the first 100-digit action evaluation.  No action, internal
gradient, finite-difference derivative, stationary classification or JSON
artifact was produced.  Therefore the run carries no positive or negative
scientific verdict about the inherited fill.

## DERIVED control defects

1. The upstream acceleration artifact is certified `10/10`; the new parser
   mistakenly required `8/8`.  Its outcome string and every other required
   upstream outcome matched.
2. State integers are encoded in layer-major order
   `code=rank+4*layer`, but one angle-lookup key used Python's default tuple
   ordering, which is rank-major.  Thus the valid triangle type
   `((0,0),(3,0),(2,1))` was stored under the differently ordered tuple
   `((0,0),(2,1),(3,0))` and raised `KeyError` before an angle was read.

Both defects are bookkeeping-only and precede the first residual.  The
geometry, action formula, dust term, precisions, finite-difference steps,
zero envelope and frozen outcome hierarchy were not exercised.

## Verdict

**DERIVED CONTROL FAILURE.**  Preserve the successful combinatorial counts,
but do not interpret them as stationarity.  Any correction must be
preregistered and limited to requiring the actual `10/10` upstream control
and using the already frozen layer-major state ordering consistently.
