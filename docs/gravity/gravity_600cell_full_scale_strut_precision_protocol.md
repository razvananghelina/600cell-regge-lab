# Protocol: high-precision resolution of the scale--strut Gram disagreement

Date: 2026-08-20

This protocol is committed before any arbitrary-precision carrier eigenvalue
is computed and before either direct LAPACK driver is compared.

## Frozen inputs

| input | SHA-256 |
|---|---|
| precision prior-art gate | `4eb4556e2c38671554db8eece1c4701fa6099dd56624c2803110a4ff9c09d015` |
| first-result note | `5753375ca2a6c4f5152f134474176501b580a1c55b7a871b3a39fa6321d82f61` |
| primary verifier | `e68105df4058f7d2ed39a6913f29e88cd9fe88e123ff52260acf698a2bd7da49` |
| first artifact | `6289b23596da28d448d1f624ecf9d9e4873ab2aa0478906dd9e90f6e13f6838d` |

The first artifact must remain `18/18` with outcome
`FULL_SCALE_STRUT_NUMERICALLY_OPEN`, exact rank 240 in both parities, zero
finite candidate mismatches, and frozen binary64 Gram discrepancies between
`0.02` and `0.04`.

## Matrix reconstruction

Reconstruct each `1560 x 240` matrix only from the frozen artifact's edge
orders, background `(lambda,rho,L0^2)` and disclosed formulas.  Do not import
or execute the primary verifier.  Require the reconstructed direct-SVD
condition number and old binary64 Gram discrepancy to reproduce the frozen
values within `1e-10` relatively.

No row or column scaling is allowed in the decisive carrier comparison.

## Known-answer calibration

Construct in binary64

```text
A_delta = diag(1,1e-7) * (1/sqrt(2)) [[1,-1],[1,1]].
```

Require:

```text
GESDD and GESVD singular values agree with (1,1e-7) within 1e-8 relative;
binary64 Gram smallest singular value has relative error > 1e-4;
80-decimal Gram singular values of the exact binary64 entries agree with
direct GESVD within 1e-8 relative.
```

If the binary64 Gram route happens not to fail on this platform, the control
fails; the carrier result may not be reinterpreted.

## Decisive carrier precision audit

For each parity:

1. compute all 240 singular values with SciPy/LAPACK `GESDD`;
2. independently compute them with SciPy/LAPACK `GESVD`;
3. convert every binary64 matrix entry to its exact integer-ratio value;
4. accumulate `G^T G` at 80 decimal digits from the sparse row supports;
5. compute all 240 symmetric eigenvalues with `mpmath.eigsy` and take their
   nonnegative square roots.

No binary64 Gram eigenvalue enters the decisive comparison.

Acceptance requires:

```text
max relative GESDD/GESVD discrepancy                < 1e-8;
max relative GESVD/80-decimal-Gram discrepancy      < 1e-8;
every 80-decimal Gram eigenvalue                     > 0;
lambda_min/(10^-80 lambda_max)                       > 1e40;
0.05 < old_binary_Gram_error/(eps*kappa_GESVD^2) < 20.
```

The large positivity margin is a precision certificate only for the frozen
binary64 matrix.  Exact rank remains certified by combinatorics, not by this
eigensolve.

## Convention and corruption controls

- Reversing the carrier row order must leave both direct spectra unchanged
  within `1e-12` relative.
- Deleting the first pole identity coefficient must reduce exact structural
  protection; the corrupted matrix's direct smallest singular value must be
  strictly smaller than the baseline value.  No claim is made that one
  deleted row forces rank loss because redundant diagonal rows may retain
  that column.
- The even and odd spectra are recorded separately; equality is not a target.

## Outcome hierarchy

1. `FULL_SCALE_STRUT_PRECISION_CONTROL_FAILED`: provenance,
   reconstruction, known-answer calibration, row-order or corruption control
   fails.
2. `FULL_SCALE_STRUT_PRECISION_RESOLVED`: both parities meet every decisive
   criterion and the old discrepancy lies on the predicted
   `epsilon*kappa^2` scale.
3. `FULL_SCALE_STRUT_PRECISION_DISAGREEMENT`: controls pass but either direct
   driver and high precision disagree, a high-precision eigenvalue is not
   positive, or its margin is insufficient.
4. `FULL_SCALE_STRUT_PRECISION_OPEN`: any remaining classified case.

Only outcome 2 removes the numerical qualifier from the accepted-background
candidate.  It does not prove the generic formula; the symbolic adversarial
gate remains mandatory.

Only this verifier and static registry guards may be run.

