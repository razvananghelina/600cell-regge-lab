# Adversarial protocol: zero corrected-strut/canonical intersection

Date: 2026-08-19

Status: post-result and target-disclosed.  The primary census is frozen in
commit `b64fd83`, artifact SHA-256
`422d8d8cb0fc0d72d842e3bf79609d4d985da6237c58e7c699b5f9cc21b65cec`.
It reports nullity zero in all 14 sectors.  This audit supplies no blind
prediction evidence; it attacks the SVD/rank calculation mechanically.

## 1. Claim under attack

For both staircase parities and every minimal `2T` sector, the matrix

```text
D = G_corrected-C_canonical
```

has full column rank `5d`.  Because the two graphs have the same literal pole
identity, this is equivalent to zero intersection of their images.

The primary route used divide-and-conquer singular values and a stacked-image
SVD.  The smallest singular value is `1.173190442e-5` in the delicate
homogeneous sector and at least `6.5697` elsewhere.  Its smallest
singular-to-uncertainty ratio is `2.79e6`.

## 2. Mechanically different decisive test

Reconstruct `D` from the same frozen geometry and action response, but do not
call the primary intersection verifier or any of its SVD/rank functions.

Compute a column-pivoted economic QR decomposition

```text
D P = Q R.
```

Solve the triangular system `R X=I` and use the rigorous algebraic inequality

```text
sigma_min(D) = sigma_min(R) >= 1/||R^-1||_F.
```

Subtract the same independently reconstructed matrix uncertainty bound from
this lower bound.  Call the result `certified_lower`.  Full column rank is
adversarially certified only if

```text
certified_lower > 100 epsilon_matrix
```

in all four derivative variants and all 14 sectors.  The decisive calculation
may use QR, triangular solves and Frobenius norms; it must not use an SVD,
eigenvalues of `D*D`, or the primary stored singular values except for a final
one-sided consistency check that the QR lower bound does not exceed the
stored singular value beyond roundoff.

## 3. Controls and convention attacks

1. Require the same literal pole identities and ranks of `G,C` before forming
   the difference.
2. Positive-intersection control `D=0` must fail the full-rank certificate.
3. Negative-intersection control given by an embedded identity must pass.
4. Apply one frozen nonsingular upper-triangular coefficient transform; the
   full-rank certificate must persist.
5. Apply a deterministic row reversal with alternating unit phases; the
   certificate must persist and its raw lower bound must agree below `2e-7`
   relative error.
6. Complex conjugation must preserve the bound below `2e-7` relative error.
7. The frozen first-diagonal source/target corruption must change at least one
   raw QR lower bound above the carrier roundoff.
8. Report the QR residual, unitarity residual, triangular condition surrogate
   and every uncertainty term.

## 4. Outcome boundary

- `CORRECTED_STRUT_ZERO_INTERSECTION_ADVERSARIALLY_CORROBORATED` only if all
  14 actual blocks have positive robust lower bounds under every variant and
  all controls pass.
- `CORRECTED_STRUT_ZERO_INTERSECTION_ADVERSARIAL_DISAGREEMENT` if a controlled
  block fails full rank.
- `CORRECTED_STRUT_ZERO_INTERSECTION_ADVERSARIALLY_OPEN` for provenance,
  reconstruction or numerical-control failures.

Corroboration closes pure-strut canonical freedom only on this fixed slab.
It does not decide mixtures with the 120 scale directions, gauge, curvature
propagation, physical tensor modes, tick duration, `c`, `G` or Planck units.

