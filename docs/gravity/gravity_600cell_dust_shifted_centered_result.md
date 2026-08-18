# Result: the shifted centered operator reproduces the `1:5` inertia

Date: 2026-08-18

## Headline

**DERIVED COMPUTATIONAL, BLIND.**  The unique centered decomposition of the
committed shifted Jacobi stencil is regular in all `56` cells and gives, for
each schedule separately,

```text
120 eigenvalues of one kinetic sign,
600 eigenvalues of the opposite sign,
0 zero-consistent,
0 open.
```

All normalized stiffness eigenvalues are real-consistent and all schedule
comparisons are robust.  The preregistered outcome is

```text
SHIFTED_CENTERED_CERTIFIED.
```

Thus the earlier finite `1:5` DeWitt-like inertia is not confined to the
first centered time slice.  This still does not identify the `120`-space as
conformal; that is the next falsifier.

## Provenance

```text
blind shifted Jacobi artifact       ee57bcc
centered protocol                   a320fe4
registered verifier                 77e4a26
```

Verifier:
`reproducible/verify_gravity_600cell_dust_shifted_centered.py`.

Artifacts:

```text
gravity_600cell_dust_shifted_centered.json
SHA-256 265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47

gravity_600cell_dust_shifted_centered.npz
SHA-256 c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8
```

Two targeted runs returned `7/7` and reproduced both artifacts byte for byte.
The full suite was not run.

## Numerical ledger

All `56` Flint determinant balls of `M_2` exclude zero.  The centered and
normalized coefficient identities hold entrywise.  The displayed midpoint
adjoint defects are small:

```text
M_2  1.23e-14 ... 2.04e-14,
N_2  6.74e-13 ... 1.19e-12,
V_2  1.70e-12 ... 3.02e-12.
```

The operational real parts span

```text
Gamma_2  -0.398712 ... -3.92e-6,
Omega_2  -7.06e-6  ...  0.319003,
```

with maximum `Omega_2` imaginary midpoint below `9.90e-14`.  All `1440`
full-multiplicity entries across the two schedules are `REAL_CONSISTENT`.

The schedule ledger is

```text
Gamma_2/Omega_2 singular spectra  14/14 SCHEDULE_ROBUST,
Omega_2 eigen spectra              7/7 SCHEDULE_ROBUST.
```

## Interpretation firewall

- **DERIVED COMPUTATIONAL:** the shifted finite centered operator exists and
  has resolved inertia `120:600` per schedule.
- **PATTERN, strengthened:** the same global `1:5` count occurs at two
  consecutive centered slices.
- **OPEN:** whether the minority space is exactly the canonical
  vertex-conformal carrier rather than another `120`-space.
- **OPEN:** whether the action-relative shape stiffness again contains the
  same `30` negative directions.
- **OPEN:** physical quotient, tensor polarizations, long-time stability,
  refinement, dispersion and speed.

The recurrence and its inertia were committed without loading any earlier
centered, conformal, shape or negative-mode result.  The next comparison may
therefore test persistence without changing this operator.

