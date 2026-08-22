# Preserved first-run failure: scalar sign-control threshold

Date: 2026-08-22.

Frozen implementation commit: `3f0dd26`.

Targeted command:

```text
/home/razvan/science/.venv/bin/python \
  reproducible/verify_gravity_600cell_finite_height_full_boundary_tangent.py
```

The first four provenance/background gates passed.  The exact scalar
generating-function control then printed

```text
det(J)                         -7 exactly
maximum good tangent error      0
good symplectic defect norm     9.66e-141
bad K_NO defect norm            4.44467
reversed-sign tangent distance  0.638876565
```

but the combined control was labelled `FAIL`.  The run was interrupted
immediately while the already accepted 43-control geometry module was being
imported.  No finite-height Regge Hessian kernel, pre-Legendre matrix,
determinant, rank, tangent or schedule comparison had been constructed.

## Cause

The frozen protocol requires the reversed pre-momentum sign to change the
scalar tangent and fail the exact expected matrix.  The implementation added
an unstated condition

```text
distance > 1
```

instead of testing resolved non-equality.  The reversed-sign matrix differs
from the exact matrix by `0.638876565`, so the protocol condition is
satisfied and the implementation-only threshold is not.

## Licensed correction

Replace only that unstated bound by

```text
distance > 1e-12.
```

This is many orders above the scalar ball/roundoff scale and tests exactly
the preregistered non-equality.  Do not change any Regge formula, derivative
step, numerical gate, outcome threshold or artifact input.  Commit the
correction before rerunning.
