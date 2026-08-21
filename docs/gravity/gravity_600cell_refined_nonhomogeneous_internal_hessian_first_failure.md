# First result: complete refined nonhomogeneous Hessian control failure

Date: 2026-08-21

Status: **10/12 CONTROL FAILED; no kernel verdict accepted.**

## Frozen provenance

- prior-art gate: `d4dc6c7`;
- target-free protocol: `fdf6f89`;
- registered implementation before first execution: `cb3e228`;
- first artifact:

  ```text
  reproducible/gravity_600cell_refined_nonhomogeneous_internal_hessian.json
  SHA-256 4a05968c68f8e6a35a1308ddf6114bb19b7106f214bfdcf798e7af2387bddec1
  ```

Only this verifier ran.  No full suite and no deferred nonlinear census ran.

## What passed

The first complete execution rebuilt all 24 labelled slabs and all actual
internal directions.  It found 96 local pentachoron patterns and 52 triangle
patterns.  All displaced stencils retained one Lorentzian Gram negative
direction, with minimum logarithm argument `0.91222028...`; the maximum local
derivative envelope was `2.21691e-31`.

Across all 24 schedules:

- every individual internal gradient passed stationarity, using at most
  `5.1811e-3` of its gate;
- raw reality and reciprocity passed, using at most `1.9652e-5` of their gate;
- the analytic product-duration tangent passed the complete local null test,
  using at most `6.6605e-10` of its gate;
- all twelve schedule/reverse pairs passed explicit layer-reversal
  congruence, using at most `1.0430e-5` of their gate.

These are retained controls, not an accepted complete-kernel result.

## The two failed controls

### Aggregate pullback

For every schedule, the maximum difference between the local orbit pullback
and the independently differentiated aggregate `10 x 10` block was exactly

```text
6.293703336268663e-10.
```

Depending on the local row-error envelope this used between `2.2504` and
`2.6369` times the preregistered gate.  The schedule-independent value is a
diagnostic pattern, not proof that the discrepancy is harmless.  Until its
entry, sign and arithmetic source are isolated, the local assembly is not
accepted.

### Corruption

The deterministic rule selected the first incidence with internal row and
column indices.  In four of twelve representatives the corresponding area
gradient was exactly zero, so the nominal `1e-4` angle-derivative corruption
changed the matrix by exactly zero.  This is a control-design bug: structural
incidence did not imply nonzero analytic coefficient.  It does not by itself
alter a scientific matrix, but the failed control cannot be retroactively
counted as passed.

## Unaccepted spectral pattern

All twelve bordered factorizations completed.  Their two eigensolves paired,
and the smallest reported absolute values lay near

```text
1.4556489532e-9.
```

Nevertheless every pair was `resolved_nonsingular = false`: the frozen
uncertainties were `8.63e-6` to `3.15e-5`, dominated by eigensolver/solve
residual terms and far larger than the observed value.  Therefore the value
is an **UNACCEPTED PATTERN**, not evidence for either a gap or an additional
null direction.

## Frozen verdict and next action

The artifact outcome is correctly

```text
LOCAL_EXTENSION_INVALID.
```

No statement about the full kernel follows.  The next admissible action is a
cheap, separately frozen diagnostic on one schedule and its `10 x 10`
pullback:

1. locate the exact discrepant entry and separate binary64 summation from a
   formula/sign error;
2. replace the ineffective corruption only after specifying a nonzero
   analytic-coefficient rule;
3. diagnose the bordered residual scaling without rerunning all twelve sparse
   factorizations.

The two-hour census will not be repeated unless this diagnostic proves that a
narrow correction can resolve the preregistered scientific question.
