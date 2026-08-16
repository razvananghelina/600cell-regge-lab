# Homothetic canonical-lapse solve: first numerical result

Date: 2026-08-16

## Provenance

- prior-art gate: `c7f3e29`;
- frozen protocol: `ded77c5`;
- implementation before evaluation: `3a90633`;
- first-run reporting failure: `4f49ebd`;
- reporting-only correction: `6d9b508`;
- corrected artifact SHA-256:
  `3ec8d39a077a595cca08cf98c17eac7b08b0e5ad01947462b63458c069bfb0f0`.

Only the targeted verifier was run.  It returns **7/7**.  The full suite was
not run.

## Verdict

Both parities return

```text
CANONICAL_LAPSE_JACOBIAN_OPEN.
```

The sole seed reproduces

```text
F = (4.318e-32, 1.616094e-9),
```

but the frozen coarse derivative calibration cannot resolve the second
Jacobian singular value, so Newton is correctly forbidden from starting.

## The recorded matrix

The operational primary matrix is identical in both schedules:

```text
[ -0.00544859760373371184   -1.69782737217922225e-8 ]
[ -582.851179621276617      -0.00136215707893697973 ].
```

Its diagnostics are

```text
singular values = (582.8511796483356..., 4.24458447860953e-9),
determinant     = -2.47396107047458e-6,
condition       = 1.37316428165e11,
epsilon         = 1.26362625715e-8.
```

The entrywise calibration passes, but the frozen rank gate requires
`s_min>100*epsilon`; here `s_min<epsilon`.  Therefore no rank-two claim is
licensed by this run.

## Interpretation

**OPEN NUMERICALLY.**  The matrix is not numerically zero-determinant, but a
nonzero determinant printed below its own calibrated error is not evidence
of invertibility.

The weak value `4.2445845e-9` coincides with the already certified weak
pseudo-constraint scale in the full canonical Legendre Jacobian
(`4.2445618e-9`).  This strongly suggests that the coarse steps, rather than
an exact rank-one identity, caused the open verdict.  That comparison is a
**PATTERN / upstream control**, not permission to reinterpret the frozen
outcome.

## Frozen next correction

A separately preregistered precision correction may change only the four
Jacobian steps to the values already successful on the same weak scale:

```text
operational: (1e-20 primary, 1e-15 shadow),
validation : (3e-20 primary, 3e-15 shadow).
```

It must retain the seed, equations, Newton rule, tolerances, branch gates and
outcome hierarchy.  If `s_min` is still inside the corrected error band, the
two-variable route remains numerically open and no solve is allowed.
