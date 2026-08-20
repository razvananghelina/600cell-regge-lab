# Preregistered correction: refined H4 stationary-fill controls

Date: 2026-08-20

First-failure record commit: `f64809b`.

This correction is frozen before the first successful action or residual
evaluation.

## Allowed changes

Exactly two semantic corrections are allowed:

1. accept the frozen projected-acceleration artifact only when
   `passed == tests == 10`, retaining its already required exact outcome;
2. canonicalize every abstract simplex and triangle state tuple by the frozen
   integer order `rank+4*layer`.  This changes only lookup-key ordering and
   must not change the enumerated vertices, simplices, triangle incidences or
   their multiplicities.

The implementation may introduce one helper for that canonical ordering and
use it at the angle-lookup construction point.

## Frozen exclusions

No input, hash, geometry formula, coordinate, action sign or branch, dust
normalization, derivative, precision, finite-difference step, error envelope,
control threshold or outcome may change.  A root, Hessian, spectrum, schedule
selection and target comparison remain forbidden.

The corrected verifier must reproduce the first run's 24 copies of
`57,600` pentachora, `149,280` triangles, 28 triangle types and zero mixed
types before reaching the variational calculation.
