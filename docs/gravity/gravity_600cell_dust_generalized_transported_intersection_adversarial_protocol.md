# Protocol: adversarial rank replication of the transported intersection

Date: 2026-08-18

Independence gate commit: `0d3a46c`.

## 1. Frozen provenance

Require exact SHA-256 hashes

```text
adversarial phase verifier
  f1cd1674af43573fd1c16b18bd37f7405093b580ce5cfa3ccad606ecb6a733cc
adversarial phase artifact
  c33615ac6d0f3133e53077f46c5ee766b9c633d4d64c32124c24839c9c84c880
earlier tangent archive
  ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d
exact intersection artifact under attack
  207cbe61bfaaf2b13d62cc3dbbb2ed5ea4931b7aab13cd47a8dd2802410c55c0
```

Replay the adversarial phase verifier and require `8/8`, outcome
`ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED`, byte-identical artifact,
32 reconstructed generalized bases and 16 tangent cells. Require the exact
intersection artifact to report `7/7` and
`TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL`; do not read its numeric singular
values.

## 2. Independent leakage rank

For each parity, sector `4,5` and all four derivative schedules, form the
canonical phase bases supplied by the adversarial path,

```text
W0 = diag(U_old,conjugate(U_old)),
W1 = diag(U_shifted,conjugate(U_shifted)),
```

load the matching earlier `T_2` midpoint and construct

```text
L = W1_perp^H T_2 W0.
```

Compute all 30 singular values with `numpy.linalg.svd`. Set the fixed numerical
floor for each cell to

```text
e = 1000 eps_machine 60 max(1,||L||_2).
```

Label every singular value only by

```text
s <= 10 e   AUDIT_SINGULAR_ZERO_CONSISTENT
s > 100 e   AUDIT_SINGULAR_NONZERO_RESOLVED
otherwise   AUDIT_SINGULAR_OPEN.
```

All 30 nonzero-resolved values certify numerical rank 30 on this independent
path. Fewer do not corroborate the exact zero intersection.

## 3. Adversarial controls

For every cell:

- positive-intersection control `T_pos=W1 W0^H`: require all 30 singular
  values zero-consistent;
- transverse control `T_neg=W1_perp W0^H`: require all 30 singular values
  nonzero-resolved;
- reverse and deterministically rephase `W0,W1`: require the complete singular
  spectrum to change by at most `10 e`;
- require finite orthonormal bases, a 30-dimensional target complement and all
  480 actual singular records.

## 4. Frozen outcome hierarchy

Use the first applicable branch:

1. provenance, replay, finiteness, dimension, census, basis-gauge or synthetic
   control fails: `ADVERSARIAL_INTERSECTION_CONTROL_FAILED`;
2. any cell has fewer than 30 nonzero-resolved singular values:
   `ADVERSARIAL_INTERSECTION_DISAGREEMENT_OPEN`;
3. all 16 cells have numerical rank 30:
   `ADVERSARIAL_INTERSECTION_ZERO_CORROBORATED`.

Only branch 3 permits the exact zero-intersection certificate to remain
accepted under project rule 4. This audit does not alter its exact error bound.

## 5. Deliverable and exclusions

Write a deterministic JSON artifact with all hashes, 16 cell records, all 480
singular records, synthetic controls, counts, outcome and status ledger.
Register before the first execution and run twice byte-identically.

Run no full suite and no unrelated verifier. Compute no graph, symplectic
restriction, propagator, dispersion, mass, inertia or limiting speed.

