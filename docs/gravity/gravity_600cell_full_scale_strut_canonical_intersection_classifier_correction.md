# Classifier correction: preserve numerical OPEN

Date: 2026-08-20

This correction is committed after the frozen `11/13` artifact and before a
corrected rerun.

## Permitted change

Keep every matrix, scaling factor, uncertainty term, singular-value driver,
threshold and corruption unchanged.  Change only the tri-state handling of
resolution-dependent controls:

- if an actual sector is resolved, the basis-change and reduced/joined
  nullities must agree or the control fails;
- if an actual sector has an `OPEN` singular value, those comparisons are
  `DEFERRED_OPEN`, not hard failures;
- the structural `C=G_strut` control passes if it resolves to nullity `5d`,
  fails if it resolves to another nullity, and is `DEFERRED_OPEN` if its
  singular spectrum contains `OPEN`;
- zero-matrix, embedded-identity, graph/pole/rank, provenance, conjugation
  and corruption controls remain hard controls.

The final hierarchy must therefore return
`FULL_SCALE_STRUT_CANONICAL_NUMERICALLY_OPEN` when all hard controls pass but
any actual or resolution-dependent control remains open.

The frozen first artifact is renamed with `_first_failure` and remains
immutable.  Only the corrected targeted verifier may be rerun; no precision
method is added by this correction and no full suite may run.

