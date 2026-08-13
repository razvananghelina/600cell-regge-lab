# Precision correction: published dust total-action differences

Date: 2026-08-13

Original protocol commit: `cc0902b`

Recorded failed run: `gravity_600cell_published_dust_control_float_failure.md`

Status: **frozen before any arbitrary-precision action evaluation**.

## Fixed mathematical problem

Do not change:

- `M_star`, `M`, `tau`, `l0`, `d` or either boundary;
- the two schedule parities;
- the Lorentzian logarithm branch;
- the gravitational or dust action;
- the 35 directions;
- relative centered-difference step `3e-6`;
- stationarity threshold `1e-7` per edge;
- the `FULL` / `POLE_ONLY` / `SOURCE_MISMATCH` labels.

## Replacement numerical control

Reimplement the **action only** for the 100 simplex orbits in `mpmath` at 60
decimal digits:

1. Cayley/edge Gram determinants and inverses;
2. signed simplex, facet and hinge volume squares;
3. the corrected complex sine/cosine angle and below-cut logarithm;
4. all orbit hinge curvatures and triangle areas;
5. the same five dust terms.

Do not import the binary64 analytic gradient into this action evaluator.
First compare its base action with the certified binary64 orbit action to
relative `3e-8`.  Then calculate all 35 centered differences at the already
frozen step.

Require:

```text
maximum difference from the analytic total orbit derivative  5e-8
maximum arbitrary-precision imaginary derivative             1e-30
maximum direct per-edge sourced residual                      1e-7.
```

Also calculate the source-free gravitational derivative and the analytic dust
derivative separately at the five poles and report their cancellation ratio.

The original binary64 finite-difference errors remain in the JSON and result
note.  They are not overwritten or relabelled as passing.

## Decision

- If all arbitrary-precision gates pass, the failed binary64 diagnostic is
  superseded as a precision limitation and the physical classification may
  be certified.
- If any arbitrary-precision gate fails, the control remains uncertified even
  if the analytic residual is small.
