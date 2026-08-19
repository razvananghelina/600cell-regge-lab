# Canonical-data carrier classifier correction

Date: 2026-08-19

## Frozen scientific question

Keep unchanged the hypotheses, candidate map, matrices, exact residual test,
negative controls, and outcome hierarchy preregistered in
`gravity_600cell_canonical_data_carrier_protocol.md`.  Keep the first failed
artifact and its SHA-256 in
`gravity_600cell_canonical_data_carrier_first_failure.md`.

## Defect

The first implementation made
`alternate_inclusion["nonzero_rows"] == 0` a construction control.  Therefore
a candidate rejected identically by the baseline and alternate exact
right-inverse graphs was reported as `CONTROL_FAILED`, preventing the already
preregistered `VERTEX_CARRIER_REFUTED` branch from being reached.

## Correction fixed before re-execution

The alternate-graph construction control is the equality of the two inclusion
decisions:

```text
(baseline nonzero rows == 0) == (alternate nonzero rows == 0).
```

The alternate candidate must still have the same derived local formulas,
exact local decomposition, and exact data rank 240.  If every construction and
negative control passes but both exact graphs reject inclusion, assign
`CANONICAL_DATA_VERTEX_CARRIER_REFUTED`.  If the two graphs disagree, assign
`CANONICAL_DATA_CARRIER_CONTROL_FAILED`.

This correction changes only verdict classification.  It must not alter a
single scientific matrix entry or residual.
