# First execution: fine-edge finite-difference control failed

Date: 2026-08-22

Status: preserved first execution, `13/14`; scientific acceptance withheld.

## Provenance

- prior-art gate: `79b612b`;
- frozen protocol: `9005af2`;
- registered implementation before execution: `077164e`;
- first artifact:

  ```text
  reproducible/gravity_600cell_nested_vertex_displacement.json
  SHA-256 dc68f3f427debb8986f2b2fa2d519a9db547f60724e26c625f7c6c3011adb50c
  ```

## What passed

Both carrier reconstructions agreed exactly.  The fine carrier contained
`2640` retained vertices and `17040` projected edge midpoints.  All midpoint
denominators were regular, all `460800` parent occurrences agreed, and the
old-vertex restriction supplied an exact index-level left inverse.  The
structural tangential rank was therefore `7920`, without an SVD.

Tangency, both determinant classes of `O(4)` covariance, all 24 schedule
controls, the target firewall and both negative controls passed.  The direct
vertex finite-difference errors were below `1.66e-10`.

These passing checks are retained as diagnostics only.  The protocol requires
all fourteen checks, so no carrier verdict is accepted from this run.

## Failed check

For the two deterministic tangent fields, the maximum fine squared-edge
derivative errors at centered steps `(2^-18,2^-20)` were

```text
(3.437332640743307e-10, 1.2806824975508668e-9),
(3.5337600901019783e-10, 1.2884276666402794e-9).
```

Every error was far below the absolute `3e-7` gate, but the errors increased
when the step was reduced and were not both below the frozen `2e-10`
roundoff exception.  The verifier correctly returned

```text
NESTED_VERTEX_PROLONGATION_DERIVATIVE_FAILED.
```

## Diagnosis and narrow admissible correction

The analytic edge derivative is formed from first-order differences, while
the control subtracts two nearly equal squared chord lengths and divides by
`2 epsilon`.  The observed approximately fourfold growth under a fourfold
step reduction is the expected binary64 cancellation signature.  This is a
diagnosis, not permission to waive the failed gate.

The narrow correction is to retain the same two steps, nonlinear normalized
map, centered formula, analytic derivative, complete edge census and frozen
tolerances, but evaluate only the finite-difference control path in
`numpy.longdouble`.  The carrier and analytic map remain binary64.  A
synthetic scalar centered-difference control must demonstrate the same
precision improvement.  If extended precision still fails the frozen
monotonicity/error rule, the derivative gate remains failed.

No action, Hessian, eigenvalue or physics target may be added during this
correction.
