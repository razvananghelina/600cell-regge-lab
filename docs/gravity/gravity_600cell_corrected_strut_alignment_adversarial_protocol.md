# Adversarial protocol: corrected strut separation from frozen dynamics

Date: 2026-08-19

Status: post-result protocol.  The primary artifact was already frozen in
commit `7ef7a7b`, with SHA-256
`5652b1371563ff11919be130af15f5b48850e2cc65a50ec35e5de85fdb587f90`.
It reports 42/42 `SEPARATED` comparisons and the outcome
`CORRECTED_STRUT_EXTREME_SELECTION_OPEN`.  Therefore this audit is not blind
and cannot add preregistration evidence.  Its purpose is to try to falsify the
numerical separation by a mechanically different subspace calculation.

## 1. Claim under attack

For each of two frozen staircase parities and seven minimal `2T` sectors, the
geometry-selected corrected 120-column strut carrier is separated from:

1. the canonical pole-Schur lift;
2. the transported largest-modulus tangent branch;
3. the transported smallest-modulus tangent branch.

The primary route used economic QR, singular values of the overlap, and
ordered complex Schur vectors for the tangent branches.  Its distances are
between `0.99779` and `0.99832`, while its largest calibrated error is below
`6.5e-8`.

The homogeneous fifth-pair modulus gap is only about `1.006 < 2`, so the
identity of the complete extreme branch remains **OPEN** by the frozen primary
protocol even if the reported candidate is numerically far from the corrected
carrier.  This audit must preserve that distinction.

## 2. Frozen inputs

Require exact hashes for:

```text
primary alignment source
d79f39380e4480aa2599d6ad0d6f56dc599268f510fe2ded1f59b9b585fb2b70

primary alignment artifact
5652b1371563ff11919be130af15f5b48850e2cc65a50ec35e5de85fdb587f90

target-blind corrected carrier artifact
e8035fb9c35ad693d1dd2adbda79485b6dd8d42bdf40a95b70a92466e47027d7

frozen old alignment source and artifact
e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047
a230a0a22c69d956b7558358d46634ad44c508326d4c34d8d7fc421aefdbcaff

full-boundary tangent numeric archive
816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b
```

The source hash was measured after freezing the primary artifact and before
writing the adversarial verifier.

## 3. Independent decisive calculation

Reconstruct the corrected carrier from its committed edge rows, and rebuild
the frozen response, lift and tangent blocks from their earlier frozen inputs.
Reusing those input operators is unavoidable because they are the objects
under audit; do not call the primary corrected-alignment verifier or any of its
QR/distance/classification functions.

For every full-rank candidate matrix `X`, construct an orthonormal basis by
the Hermitian polar route

```text
G = X* X = V diag(mu) V*
Q = X V diag(mu^(-1/2)) V*.
```

Construct `P = Q Q*`.  Obtain the projector distance from the largest absolute
eigenvalue of the explicitly Hermitian matrix `P_corrected - P_target`.  Do
not use QR or an SVD of an overlap matrix for this decisive value.

For the two tangent targets, diagonalize the frozen tangent matrix directly,
sort eigenvectors by eigenvalue modulus, take the fixed `5d` largest or
smallest vectors, and transport them through the response.  Do not use ordered
Schur vectors.  Retain and report the frozen modulus gaps.

## 4. Fixed tests and thresholds

Evaluate exactly the same 42 named comparisons.  Require:

- all polar Gram eigenvalues positive and all orthonormality residuals below
  `2e-11`;
- every adversarial projector distance greater than `0.99`;
- every distance differs from the corresponding primary distance by less
  than `2e-6`;
- the last one-dimensional sector has both modulus gaps below `2`, while all
  other sectors retain gaps above `2`;
- even/odd distances differ by less than `2e-7` without combining them.

The `0.99` lower bound is deliberately far from both zero and the maximum
reported numerical uncertainty.  Because it was chosen after seeing the
primary result, it is a robustness threshold only, not evidence of a blind
prediction.

## 5. Known controls and convention attacks

1. **Positive control:** multiply every corrected basis by a frozen
   deterministic nonsingular upper-triangular matrix.  The polar projector
   distance to the original span must be below `2e-11`.
2. **Negative control:** add `0.05` times a deterministic vector in the
   orthogonal complement to one column.  The span must change by more than
   `1e-3`.
3. Complex-conjugate every sector basis and all matrices consistently; the
   42 distances must be invariant below `2e-7`.
4. Reverse the source/target coefficients on the frozen first diagonal; at
   least one corrected-sector projector and one target distance must change
   above `1e-10`.
5. Report condition numbers, orthonormality residuals and direct-eigenvector
   residuals.  Any failure leaves the result **OPEN**.

## 6. Outcome boundary

- `CORRECTED_STRUT_SEPARATION_ADVERSARIALLY_CORROBORATED` only if all controls
  pass and all 42 mechanically different distances satisfy the thresholds.
- `CORRECTED_STRUT_SEPARATION_ADVERSARIAL_DISAGREEMENT` if a controlled
  calculation finds an identification or a distance disagreement.
- `CORRECTED_STRUT_SEPARATION_ADVERSARIALLY_OPEN` for unresolved numerical,
  conditioning, provenance or convention failures.

Even corroboration establishes only a negative subspace statement for these
three frozen targets.  It does not identify gauge, curvature propagation,
gravitons, a physical instability, a tick, `c`, `G` or a Planck scale.
