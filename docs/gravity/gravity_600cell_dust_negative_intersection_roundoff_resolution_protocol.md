# Protocol: fixed-input resolution of the negative-intersection audit

Date: 2026-08-18

Framing gate commits: `cf64f34`, `ba42d40`.

Status: **TARGET-DISCLOSED, PREREGISTERED BEFORE ANY FIXED-INPUT
HIGH-PRECISION LEAKAGE SPECTRUM.**

## 1. Frozen provenance

Require these exact SHA-256 hashes:

```text
binary negative-fiber source
  f462e507500d7f02ecf799f0d4b320e05795216a36a0d10eb908d6dc67b48181
binary negative-fiber artifact
  d630bf07066f88c35eee5a62a80ec1f43399a95ea882a43528289220c67f4599
binary tangent archive
  ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d
primary exact-intersection artifact
  c490431bdaeae3026692cd358f60d0b47ef5d63aa59217e400daac807ed21be0
first adversarial source
  6aa7e841d31bdc87568a6e4370ed334b6f5c09884669ee13fcd68d46ea4b3162
first adversarial artifact
  d5074507326bb981ad7573bd562c1aa9f0af4e1eb6b6924e3ac959a5fa1d3340
```

Replay the binary negative verifier and require `8/8` with its accepted
outcome and unchanged artifact.  Require the primary artifact to report
`10/10` and `NEGATIVE_TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL`.  Require
the first audit artifact to report `11/11` and
`ADVERSARIAL_NEGATIVE_INTERSECTION_DISAGREEMENT_OPEN`.  Do not import the
primary verifier or read its numeric singular values.

## 2. Exact fixed binary inputs

For every real and imaginary float component in the binary projectors and
tangent midpoint, use `float.as_integer_ratio()` and form the exact dyadic
rational in `mpmath`.  No decimal serialization or source-ball radius enters
this bridge.

Run the entire construction independently at decimal precisions `100` and
`140`.  At each precision use

```text
arithmetic_floor = 1e-75.
```

This retains at least `25` guard digits at the lower precision and is frozen
independently of any observed singular value.

## 3. Projector spectral ranges and residual bound

For every binary projector midpoint `P`, form `(P+P*)/2`, diagonalize it with
`mpmath.eighe`, and let `U_-` and `U_+` contain the bottom and top `15`
eigenvectors.  Require:

- `lambda_14 < 1e-10` and `lambda_15 > 1-1e-10`;
- `||U*U-I||_2 < 1e-70` at both precisions;
- finite ordered eigenvalues and positive gap
  `g=lambda_15-lambda_14`;
- top eigenspace residual
  `r=||P U_+ - U_+ diag(lambda_15,...,lambda_29)||_2`.

Define only when `g>2r`

```text
eta = 2r/(g-2r) + arithmetic_floor.
```

Build the phase fiber and its explicit spectral complement as

```text
W     = diag(U_+,conjugate(U_+)),
W_perp= diag(U_-,conjugate(U_-)).
```

## 4. Square leakage and complete error

For old/shifted errors `eta_0,eta_1`, exact dyadic tangent midpoint `T` and

```text
L = W_1,perp* T W_0,
```

use the fixed-input perturbation bound

```text
e_L = (eta_0+eta_1+eta_0 eta_1)||T||_2
    + arithmetic_floor max(1,||T||_2).
```

This controls only invariant-subspace reconstruction and high-precision
arithmetic relative to the exact stored binary inputs.

Compute all `30` singular values three ways:

1. direct `mpmath.svd`;
2. square roots of the ordered eigenvalues of `(L*L+(L*L)*)/2`;
3. for the smallest value, `1/||L^-1||_2`.

Require the complete SVD/Gram spectra and the inverse-norm smallest value to
agree within `10 e_L`.  Label each direct singular value only by

```text
s <= 10 e_L   FIXED_SINGULAR_ZERO_CONSISTENT
s > 100 e_L   FIXED_SINGULAR_NONZERO_RESOLVED
otherwise     FIXED_SINGULAR_OPEN.
```

At both precisions require all `30` values nonzero-resolved.  Compare the
complete `100`- and `140`-digit spectra and require each difference at most
`10(e_100+e_140)`.

## 5. Controls and census

At each precision and in every cell require:

- `T_full=W_1 W_0*` gives all `30` leakage values zero-consistent;
- `T_zero=W_1,perp W_0*` gives all `30` leakage values
  nonzero-resolved;
- deterministic reversal/rephasing of source fiber and target complement
  preserves the complete spectrum within `10e_L`;
- all bases are `60 x 30`, leakage is `30 x 30`, values are finite and sorted.

Require exactly `32` projector records per precision, `16` cell records, and
`960` singular records (`16 x 30 x 2` precisions).

## 6. Frozen outcome hierarchy

Use the first applicable branch:

1. provenance, replay, dyadic conversion, projector split/residual,
   orthonormality, SVD/Gram/inverse agreement, precision stability, controls,
   finiteness or census fails:
   `NEGATIVE_INTERSECTION_ROUNDOFF_RESOLUTION_CONTROL_FAILED`;
2. any cell/precision has fewer than `30` nonzero-resolved values:
   `NEGATIVE_INTERSECTION_ROUNDOFF_DISAGREEMENT_REMAINS_OPEN`;
3. all `16` cells have rank `30` at both precisions:
   `NEGATIVE_INTERSECTION_ROUNDOFF_DISAGREEMENT_RESOLVED`.

Branch 3 says the first audit's open classification arose from its global
float64 envelope, not a rank loss in its fixed binary data.  Together with
that audit's independent spectrum and the separate source-certified primary
result, it permits consolidation under project rule 4.  It does not alter the
first audit artifact or claim that binary midpoints enclose the physical
source.

## 7. Deliverable and exclusions

Write a deterministic JSON artifact with all hashes, per-precision projector
and cell controls, all `960` values, precision differences, outcomes and a
status ledger.  Register before first execution and run twice byte-identical.

Run no full suite and no unrelated verifier.  Compute no graph, Lagrangian
restriction, propagator, dispersion, graviton, mass, inertia or limiting
speed.
