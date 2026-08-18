# Protocol: adversarial replication of the negative-intersection rank

Date: 2026-08-18

Independence gate commit: `8d532e8`.

Status: **PRIMARY TARGET DISCLOSED, PREREGISTERED BEFORE ANY AUDIT SPECTRUM.**

## 1. Frozen provenance

Require exact SHA-256 hashes:

```text
binary negative-fiber source
  f462e507500d7f02ecf799f0d4b320e05795216a36a0d10eb908d6dc67b48181
binary negative-fiber artifact
  d630bf07066f88c35eee5a62a80ec1f43399a95ea882a43528289220c67f4599
earlier tangent midpoint archive
  ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d
primary exact-intersection artifact under attack
  c490431bdaeae3026692cd358f60d0b47ef5d63aa59217e400daac807ed21be0
```

Replay the binary negative-fiber verifier and require `8/8`, outcome
`NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED`, byte-identical artifact and `32`
binary projectors.  Require the primary artifact to report `10/10` and
`NEGATIVE_TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL`, but do not read or
compare its numeric singular values.

## 2. Independent phase bases and leakage

For each old/shifted binary projector `P`, diagonalize `(P+P*)/2` anew with
`scipy.linalg.eigh` and use the top `15` eigenvectors as `U`.  Require the
bottom `15` projector eigenvalues below `1e-10`, the top `15` above
`1-1e-10`, and `||U*U-I||_2<1e-12`.

Build

```text
W = diag(U,conjugate(U))
```

and require `||W*W-I||_2<1e-12`.  Load the matching earlier binary tangent
midpoint `T_2`.  Obtain the target complement from complete QR,

```text
[W_1,W_1_perp] = qr(W_1, mode=complete),
L = W_1_perp* T_2 W_0.
```

Compute all `30` singular values with `numpy.linalg.svd`.  Freeze the cell
roundoff floor

```text
e = 1000 eps_machine 60 max(1,||T_2||_2,||L||_2).
```

Assign only

```text
s <= 10 e   AUDIT_SINGULAR_ZERO_CONSISTENT
s > 100 e   AUDIT_SINGULAR_NONZERO_RESOLVED
otherwise   AUDIT_SINGULAR_OPEN.
```

All `30` values nonzero-resolved means independent numerical rank `30`.
Fewer do not corroborate the primary certificate.

## 3. Adversarial controls and stresses

For every cell require:

1. `T_full=W_1 W_0*` gives all `30` leakage singular values
   zero-consistent;
2. `T_zero=W_1_perp W_0*` gives all `30` values nonzero-resolved;
3. reversing columns and multiplying both phase bases by deterministic unit
   phases changes the actual spectrum by at most `10e`;
4. replacing the complete-QR complement by
   `scipy.linalg.null_space(W_1*)` changes the spectrum by at most `10e`;
5. replacing NumPy's default SVD by SciPy `gesvd` changes the spectrum by at
   most `10e`;
6. swapping the `(q,p)` blocks simultaneously in source, target and tangent
   changes the spectrum by at most `10e`.

Require finite `60 x 30` phase bases, `60 x 30` complements, `30 x 30`
leakage matrices, sorted spectra and exactly `480` actual singular records.

## 4. Frozen outcome hierarchy

Use the first applicable branch:

1. provenance, replay, projector split, orthonormality, dimensions,
   finiteness, census, synthetic controls or any stress fails:
   `ADVERSARIAL_NEGATIVE_INTERSECTION_CONTROL_FAILED`;
2. any actual cell has fewer than `30` nonzero-resolved values:
   `ADVERSARIAL_NEGATIVE_INTERSECTION_DISAGREEMENT_OPEN`;
3. all `16` cells have numerical rank `30`:
   `ADVERSARIAL_NEGATIVE_INTERSECTION_ZERO_CORROBORATED`.

Only branch 3 permits consolidation of the exact primary result under project
rule 4.  This binary64 audit does not replace or strengthen the exact error
bound.

## 5. Deliverable and exclusions

Write a deterministic JSON artifact containing all hashes, basis controls,
`16` cell records, all `480` actual singular records, synthetic controls,
stress differences, counts, outcome and status ledger.  Register the verifier
before first execution and run twice with a byte-identical artifact.

Run no full suite and no unrelated verifier.  Compute no graph, Lagrangian
restriction, propagator, dispersion, graviton, mass, inertia or limiting
speed.
