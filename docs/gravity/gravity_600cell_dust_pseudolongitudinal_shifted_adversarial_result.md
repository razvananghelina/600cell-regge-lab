# Result: direct shifted residuals persist, but the preregistered rank conjunction is OPEN

Date: 2026-08-19

## Literal preregistered outcome

The direct adversarial verifier returned

```text
SHIFTED_PSEUDOLONGITUDINAL_DIRECT_OPEN
18/18 checks passed
```

This outcome is preserved without reinterpretation.  Its artifact is commit
`f0f0625` with SHA-256

```text
9e9f7253fd10422f3534914fae020857162862123fd4eae889e3570083552179.
```

The preregistered protocol is commit `1bfe9e9`; the registered verifier is
commit `3ea94a0`.

## What the independent construction found

The adversarial route did not load the shifted centered `M,V` archive.  It
rebuilt slabs 2 and 3 from local 4-simplex Hessians at 100 decimal digits,
reconstructed exact golden-ratio vertices and all 720 edges, and used SVD
rather than the primary QR path.  It opened the primary residual artifact
only after its 16-cell direct census.

Every direct cell retained

```text
stiffness inertia                  15 negative + 10 positive
rho_span                           NONZERO_RESOLVED
rho_comm                           NONZERO_RESOLVED
exact residual norm inequalities  PASS
```

For the operational-primary cells, the direct values are approximately

```text
rho_span = 0.077762667706
rho_comm = 0.009018893331
```

and the direct/archived ratios over all 16 cells lie in

```text
rho_span  1.000000021252 ... 1.000000021258
rho_comm  1.000000021273 ... 1.000000021279.
```

Thus the two independent routes agree on the nonzero residuals to about two
parts in `10^8`.  This agreement is **DERIVED COMPUTATIONAL**.  It is not a
claim that the residual is conserved.

## Why the frozen outcome remained OPEN

The protocol conjunctively required an auxiliary augmented-rank classifier
to return a rank strictly greater than 15.  It returned 15 in all 16 cells:

```text
sixteenth augmented singular value  about 2.96331418e-5
frozen augmented threshold           about 4.25805040e-5
```

The threshold therefore did not resolve the additional rank, even though the
span residual itself was resolved nonzero.

This is a classifier mismatch.  For `X=B L`, `Y=A L` and the orthogonal
projector `P_X` onto `im X`, exact linear algebra gives

```text
rank([X,Y]) = rank(X) + rank((I-P_X)Y).
```

Consequently, once `rank(X)=15` is resolved, a nonzero
`||(I-P_X)Y||` and augmented rank greater than 15 are not independent
requirements.  The first implies the second.  The frozen rank threshold
instead multiplies a dimensionless residual error by the largest singular
value of `[X,Y]`, whose scale is dominated by `X`; it is not a calibrated
error bound for the small orthogonal block.  In the direct protocol the
formula also adds an explicit binary rounding floor to `matrix_error`, which
already contains the same binary rounding term.  The resulting local floor
is almost exactly twice that of the primary audit.

The OPEN outcome is therefore an honest negative about the preregistered
conjunction, not evidence that the direct residual vanished.  No threshold is
changed in the frozen verifier.

## Status ledger

- **DERIVED COMPUTATIONAL:** two mechanically different constructions find
  both shifted pseudo-longitudinal residuals nonzero in all 16 cells.
- **DERIVED / STRUCTURAL:** exact nonzero span leakage is equivalent to an
  augmented image larger than `im(BL)`; the two tests are not independent.
- **OPEN UNDER THE FROZEN HIERARCHY:** the direct audit's auxiliary numerical
  rank threshold did not resolve rank greater than 15.
- **PATTERN:** current and shifted normalized residuals are nearly equal.
  This can partly reflect scale invariance and a near-homothetic recurrence;
  two samples do not establish conservation.
- **OPEN:** curvature/refinement scaling, continuum gauge recovery, physical
  instability, propagation, polarization and speed.

## Required correction

Any correction must be committed separately and target-disclosed.  It may
use the exact rank identity and the already frozen, independently rebuilt
span residual, but it may not overwrite the OPEN artifact or silently lower
the augmented singular-value threshold.

