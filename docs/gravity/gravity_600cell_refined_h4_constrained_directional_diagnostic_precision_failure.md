# Directional diagnostic control failure: directions parsed at default precision

Date: 2026-08-21

After the exact complex-parser repair, the diagnostic completed with

```text
14/15 PASS
REFINED_H4_DIRECTIONAL_DIAGNOSTIC_CONTROL_FAILED
```

and artifact

```text
reproducible/gravity_600cell_refined_h4_constrained_directional_diagnostic.json
SHA-256 4d4ec1489bdf2119750674ae0cb1caf724ea88f86245dd51398a5796b329a953.
```

The load-bearing reproduction control failed: the first Richardson estimate
differed from the frozen primary value by as much as `1.145e-9`.  Inspection
of the code shows that the stored 65--70 digit basis, lifts and quadratic
values were converted with `mp.mpf` before entering either `mp.workdps`
context.  They were therefore created at the process default precision
(approximately 15 decimal digits) and only subsequently carried into the
140/180-digit action evaluations.

This explains why the apparent Richardson error ratios were approximately one
and why none of the extrapolated values matched.  Those downstream numbers
are invalid as a scientific diagnostic because their input direction was not
the frozen high-precision direction.  The polynomial control passed because
it did not use that reconstruction path.

The failure is preserved.  No truncation or Hessian/action verdict is
accepted from this artifact.

