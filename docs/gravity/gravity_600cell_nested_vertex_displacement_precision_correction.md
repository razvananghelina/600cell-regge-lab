# Precision correction: fine-edge derivative control

Date: 2026-08-22

Status: frozen after the preserved `13/14` first execution and before changing
or rerunning the verifier.

## Frozen failed inputs

```text
registered first implementation
  reproducible/verify_gravity_600cell_nested_vertex_displacement.py
  SHA-256 d60237118b4c8ba4452100dd77b4462b092ff4b5210adc095a36d8a694a40f80

preserved first artifact
  reproducible/gravity_600cell_nested_vertex_displacement.json
  SHA-256 dc68f3f427debb8986f2b2fa2d519a9db547f60724e26c625f7c6c3011adb50c
```

The first execution and its scientific non-acceptance are recorded in
`gravity_600cell_nested_vertex_displacement_first_failure.md`.

## Allowed correction

Change only the centered finite-difference control path to
`numpy.longdouble`:

1. cast the coarse positions, tangent probes and epsilon to `longdouble`;
2. normalize rows by an explicit `sqrt(sum(x*x))` expression;
3. construct the perturbed fine positions from the same singleton/edge key
   arrays and the same normalized midpoint formula;
4. form the same centered vertex and squared-edge derivatives;
5. compare them with the unchanged analytic binary64 derivatives after
   casting those analytic values to `longdouble`.

Do not change:

- either carrier or key order;
- either deterministic tangent field;
- the two steps `(2^-18,2^-20)`;
- the centered formula;
- any tolerance or monotonicity rule;
- the analytic prolongation;
- the fourteen-check count or outcome hierarchy.

Add inside the existing edge-derivative check a scalar cancellation control

```text
f(epsilon) = (1 + epsilon)^2,
f'(0) = 2,
```

evaluated by the same centered formula in binary64 and long double at
`epsilon=2^-20`.  Require the long-double absolute error to be no greater
than the binary64 error.  This does not add a fifteenth top-level check.

## Acceptance boundary

The corrected run advances only if all original fourteen checks pass with
the frozen thresholds and monotonicity exception.  Run it twice and require
byte-identical artifacts.  Otherwise retain
`NESTED_VERTEX_PROLONGATION_DERIVATIVE_FAILED`.

No full suite, action, Hessian, spectral target or physics interpretation is
authorized by this correction.
