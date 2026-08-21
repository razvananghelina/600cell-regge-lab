# First execution failure: synthetic block used the production dimensions

Date: 2026-08-21

The registered verifier in commit `69ace62` stopped before constructing any
new 600-cell Hessian and before writing an artifact.  The frozen provenance,
upstream, topology, on-shell, Lorentzian-branch and algebraic-basis checks had
passed.

The synthetic `2+2` control called `restricted_response`, whose block splitter
was deliberately specialized to the production `12+10` dimensions.  It
therefore raised

```text
IndexError: matrix index out of range
```

inside `hessian_blocks`.  This is an implementation error in the independent
small control, not a scientific outcome.  No constrained response matrix,
schedule class, reversal comparison or directional action check was computed.

The admissible correction is restricted to evaluating the synthetic example
directly from its already frozen `A,B,C,p,q` blocks.  The expected scalar
`18`, incompatible residual `3`, production block splitter, physical inputs,
bases, thresholds, outcomes and all other code must remain unchanged.

