# Protocol: adversarial replication of generalized phase transport

Date: 2026-08-18

Prior-art and independence gate commit: `f071366`.

## 1. Frozen inputs

Require these exact files and SHA-256 hashes:

```text
direct verifier
  01479fcaa7e5354ea3bb72306ac8cd433a87b11a539f912075d69273a014b510
direct artifact
  8ded406366dbf291da02dfbf995c4e37036cc6ce745d9240d14905664ba6042a
two-step tangent verifier
  c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717
two-step tangent metadata
  f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc
two-step tangent archive
  ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d
accepted phase-transport artifact under attack
  45eb9a3e80ead758d9b3c2f8e1eccff44b06e2759251ab00c447aa53e6705743
```

Require the direct artifact to report `13/13` and
`DIRECT_GENERALIZED_COMMON_BUNDLE_RESOLVED`; replay the direct verifier and
require the artifact to remain byte-identical. Require the tangent metadata
to report `16/16`, `TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED`, and the archive
hash above. Require the result under attack to report `7/7` and
`GENERALIZED_PHASE_TRANSPORT_REFUTED`; do not read its numeric residuals.

## 2. Independent generalized bases

For every old/shifted time, parity, sector `4,5` and all four schedules:

1. reconstruct the five-dimensional conformal carrier from the raw incidence
   SVD;
2. compute the null space of `U^H M` with `scipy.linalg.null_space`;
3. restrict `B=-M` and `A=-V` to that 25-dimensional shape carrier;
4. factor `B=L L^H` explicitly with Cholesky;
5. form `C=L^-1 A L^-H` by triangular solves and diagonalize its Hermitian
   symmetrization with `numpy.linalg.eigh`;
6. lift the 15 negative eigenvectors and QR-orthonormalize them.

Require exactly `15` negative and `10` positive eigenvalues, positive `B`,
orthonormality defect at most `1e-10`, generalized eigen-residual at most
`1e-9`, and distance at most `1e-8` from the earlier direct-precision
projector as a construction control. The earlier projector is not used in
the transport test.

## 3. Basis/image test

Load, but do not reconstruct, each committed second-slab tangent midpoint
`T_2`. Require all 16 midpoints finite, `60 x 60`, with finite archive-radius
Frobenius norm. For each cell make the canonical phase bases

```text
W0 = diag(U_old, conjugate(U_old)),
W1 = diag(U_shifted, conjugate(U_shifted)).
```

QR-orthonormalize `Y=T_2 W0` to obtain `Z`. Require numerical image rank 30
using the fixed relative threshold `1e-12`. Compute without the accepted
verifier's four-block formula:

```text
sine_max = || W1_perp^H Z ||_2,
relative_least_squares = ||Y-W1(W1^H Y)||_2 / ||Y||_2,
projector_crosscheck = ||Z Z^H-W1 W1^H||_2.
```

The first value is decisive. The projector expression is only a redundant
crosscheck and must agree with `sine_max` to `1e-10`.

Frozen labels:

```text
sine_max <= 1e-10   INDEPENDENT_CLOSED
sine_max >  1e-6    INDEPENDENT_NONCLOSING
otherwise           INDEPENDENT_OPEN
```

## 4. Positive, negative and convention controls

For every cell:

- positive map `T_pos=W1 W0^H`: require `sine_max <= 1e-10`;
- negative map: replace the first column of `W1` by the first orthogonal-
  complement column and map `W0` to that basis; require
  `sine_max >= 0.5`;
- reverse and deterministically rephase source and target basis columns;
  require the decisive value to change by at most `1e-10`.

Also report, without affecting the outcome, the same calculation for
`diag(U,U)` and `diag(conjugate(U),U)`. These are convention stress tests,
not alternative definitions of the frozen cotangent lift.

## 5. Frozen outcome hierarchy

Use the first applicable branch:

1. provenance, replay, basis, rank, finiteness, crosscheck or synthetic control
   failure: `ADVERSARIAL_PHASE_TRANSPORT_CONTROL_FAILED`;
2. any canonical cell is `INDEPENDENT_CLOSED`:
   `ADVERSARIAL_PHASE_TRANSPORT_CONTRADICTED`;
3. no cell is closed but any is open:
   `ADVERSARIAL_PHASE_TRANSPORT_DISAGREEMENT_OPEN`;
4. all 16 canonical cells are nonclosing:
   `ADVERSARIAL_PHASE_TRANSPORT_REFUTATION_CORROBORATED`.

Only branch 4 permits the high-precision refutation to remain accepted under
project rule 4. This audit is **STRUCTURAL INDEPENDENT CORROBORATION**, not a
second exact proof.

## 6. Deliverable and exclusions

Write a deterministic JSON artifact containing provenance, 32 basis controls,
16 decisive records, 32 convention-stress records, all synthetic controls,
counts and outcome. Register the verifier before its first execution and run
it twice with a byte-identical artifact.

Run no full suite and no unrelated verifier. Do not compute transported
intersections, graphs, propagators, roots, dispersion, masses or speeds in
this audit.

