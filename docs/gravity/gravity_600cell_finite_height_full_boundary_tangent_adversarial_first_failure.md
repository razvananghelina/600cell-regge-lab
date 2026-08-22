# First adversarial finite-height full-boundary tangent failure

Date: 2026-08-22

## Frozen provenance

- adversarial protocol commit: `3be6eb6`
- registry commit: `b3ab177`
- adversarial implementation commit: `fb80285`
- targeted verifier:
  `reproducible/verify_gravity_600cell_finite_height_full_boundary_tangent_adversarial.py`
- the full suite was not run

## What happened

The first targeted scientific run completed both independently assembled dense
real Hessian routes.  It then classified the even full pre-Legendre matrix as
regular and passed the even hostile omission of the actual `K_NO` term.  The
odd solve also completed, but the script stopped before assigning its label:

```text
TypeError: '<' not supported between instances of 'float' and 'NoneType'
```

The hostile `K_NO` omission was deliberately evaluated only on the even
schedule.  The common bookkeeping expression nevertheless tried to compare
the absent odd hostile defect (`None`) with its numerical envelope.

Elapsed wall time was 1:43.31 and maximum resident memory was 1,032,904 kB.

## Scientific status

- **DERIVED COMPUTATIONAL, PARTIAL:** both dense 2280 by 2280 Hessian
  assemblies were finite, real and inside the preregistered raw reciprocity
  envelope.
- **DERIVED COMPUTATIONAL, PARTIAL:** the even 1560 by 1560 pre-Legendre
  matrix was classified regular, with normalized smallest singular value
  `4.0880740195610725e-4` against gate `3.651e-7`.
- **DERIVED COMPUTATIONAL CONTROL, PARTIAL:** omitting the even `K_NO` term
  produced symplectic defect `5.919e1` against gate `5.366e-1`.
- **OPEN:** the odd rank label, both final symplectic labels, schedule
  robustness and comparison with the primary representation-theoretic route.

No scientific threshold, derivative step, normalization or outcome hierarchy
is changed.  The only licensed correction is to mark the hostile control as
not applicable on the odd route instead of comparing `None` numerically.

