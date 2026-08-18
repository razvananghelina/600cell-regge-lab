# Stationary-root enumeration: test-counter reporting failure

Date: 2026-08-16

## Provenance

- prior-art gate: `eecc80e`;
- frozen protocol: `07083cc`;
- implementation before evaluation: `0c48e22`;
- first artifact SHA-256:
  `da0897033e4a0203a20b67d9fc051c6f33406b93f3f82c87681c2411bcc91f6a`.

The targeted verifier evaluated the frozen grid and printed five `PASS` lines,
but ended with `RESULT: 3/5` and exit code 1.

## Cause

At module scope, the sign-bracket refinement assigned its local boolean to the
name `passed`.  That name is also the global test counter.  The assignment
replaced the counter with `True` before the final two checks, producing
`1+2=3` despite all five checks passing.

This is a reporting/harness bug.  It does not alter the grid, candidates,
bisections, equations, roots, momenta, derivative matrices or outcome.

## Frozen scientific payload

The artifact records outcome `STATIONARY_ROOTS_ENUMERATED`, one sign bracket,
one near-zero node and two complete roots.  These values remain provisional
until the identical run is reproduced with a correct counter.

The sole permitted patch renames the refinement-local boolean.  No scientific
constant or criterion may change.  Only the targeted verifier is rerun.
