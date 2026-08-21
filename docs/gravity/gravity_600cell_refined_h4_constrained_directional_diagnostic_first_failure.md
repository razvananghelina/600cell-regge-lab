# First directional-diagnostic execution failure: stored complex parser

Date: 2026-08-21

The registered diagnostic in commit `42b3b88` completed all twelve action
ladders at both 140 and 180 digits, but stopped before evaluating the frozen
outcome or writing an artifact.  Provenance, upstream, definition, topology,
direction reconstruction and displacement controls had passed.

The optional reproduction control attempted to parse a stored `mpmath`
complex string such as

```text
(-329217... - 2.35e-116j)
```

by calling `mp.mpc(text)`.  That constructor does not accept the serialized
parenthesized `a +/- bj` form and raised `ValueError`.

No truncation, extrapolated-match or scientific verdict was assigned.  The
only admissible correction is an exact decimal parser that splits the already
stored real and imaginary fields and passes them separately to `mp.mpf`.
Steps, precisions, action evaluations, extrapolation, thresholds, controls and
outcomes must remain unchanged.

