# Preserved first run failure: branch-history serialization precision

Date: 2026-08-22

Implementation commit: `fa98abd`.

Status: **CONTROL FAILURE PRESERVED; NO SECOND-SLAB TANGENT RESULT ACCEPTED.**

## Failure

The first targeted scientific execution passed exact provenance, registry,
accepted-input, symbolic scale, geometry and carrier controls.  It then
reported

```text
[FAIL] deterministic bisection reconstructs the committed branch-B two-slab history
q1=9.618002653341898097389452521
q2=31.2792236208252636198289982
committed_error=1.67408e-59
```

The failure was diagnosed before any pre-Legendre matrix or tangent was
completed.  The process was deliberately interrupted because the frozen
outcome would already have been `SECOND_FULL_BOUNDARY_TANGENT_CONTROL_FAILED`.
No JSON or NPZ result artifact was produced.

The run had begun the even normalized and physical local derivatives.  Their
printed `SCALE_LIFT_CONFIRMED` label is discarded as an incomplete preliminary
value and is not evidence for the scientific outcome.

## Exact cause

The primary composition verifier serializes values with

```python
mp.nstr(value, 60)
```

and its accepted artifact therefore contains approximately 60 significant
digits.  The history-provenance correction incorrectly required a numerical
difference below `1e-65` against that 60-digit string.  The observed
`1.67408e-59` discrepancy is the unavoidable serialization remainder, not a
branch or equation disagreement.

The mechanically different composition artifact stores a 180-digit run and
does support the original `<1e-65` numerical comparison.

## Frozen correction before rerun

The rerun must use both controls at their actual information content:

1. reconstruct `q2,h2,r2` without using either artifact as a seed;
2. require exact equality of `mp.nstr(reconstructed_value, 60)` with the three
   primary strings;
3. retain the `<1e-65` numerical comparison against the 180-digit
   adversarial branch-B record and against the accepted first-tangent
   background;
4. retain every residual, bracket, positivity, junction, derivative,
   scale-lift and outcome gate unchanged.

This does not loosen a scientific tolerance.  It replaces an impossible
comparison with an exact test of all digits that the primary artifact stores,
while preserving the stronger numerical test on the independent high-
precision witness.
