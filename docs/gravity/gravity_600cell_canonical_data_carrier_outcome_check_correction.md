# Canonical-data carrier outcome-check correction

Date: 2026-08-19

## Frozen scientific result

Keep unchanged every hypothesis, matrix, candidate map, exact residual,
negative control, and the baseline/alternate decision correction.  The exact
vertex scale/lapse candidate is rejected on 3600 face rows under both
right-inverse graphs.

## Remaining harness defect found by source audit

The verifier calls `check(..., inclusion_ok)` on the target proposition itself.
Consequently, the preregistered and scientifically valid
`CANONICAL_DATA_VERTEX_CARRIER_REFUTED` outcome necessarily leaves one red
test and exits nonzero.  A falsification target is not a construction control;
either truth value must be accepted if it was computed completely and the
outcome hierarchy records it.

## Correction fixed before the next execution

Replace the target-valued check by an inclusion-scan accounting check.  For
each representative, require

```text
0 <= nonzero_rows <= number_of_complete_rows
(first_nonzero is None) == (nonzero_rows == 0)
```

for both the baseline and alternate graphs.  Keep `inclusion_ok` unchanged as
the scientific decision fed to the outcome hierarchy.  Add the accounting
condition to `controls_ok`.

Thus:

- complete agreeing zero residuals may reach the positive branch;
- complete agreeing nonzero residuals reach `VERTEX_CARRIER_REFUTED`;
- incomplete accounting or disagreement remains `CONTROL_FAILED`.

This is a verdict-harness repair only.  It must not change any scientific
entry in either preserved failed artifact.
