# Adversarial protocol: geometric replication of the complete Coxeter census

Date: 2026-08-21

Status: frozen after the primary all-schedule result `d0adbac` and before
constructing any adversarial block for a newly tested schedule.

## 1. Claim under attack

The primary census reports, under the fixed refined-slab hypotheses,

```text
ker(C_s) = span(n_s)
```

for all twelve time-reversal representatives.  Schedule 0 was already
adversarially corroborated in `66e47a7`.  The exact claim under attack here is
the same statement for representatives

```text
1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 14.
```

The primary spectra and their repeated weakest `k=1` values are known before
this protocol.  They are targets for replication, not discovery evidence.
No tolerance, sector, schedule or method may be selected from their values.

## 2. Independence boundary

This replication attacks the group action, Fourier reduction and numerical
spectrum by a mechanically different route.  It does **not** independently
derive the Regge Hessian itself: both routes use the same frozen local action,
coordinates, masses and binary64 source-matrix assembly.  Therefore a
positive result corroborates the kernel of those exact frozen matrices; it
does not exclude a common error in the upstream physical action.

The primary coloured-graph automorphism and entrywise Fourier accumulation
are forbidden.  The adversarial construction must instead use:

1. ambient `R4` reflection matrices derived geometrically from the normalized
   centres of a base barycentric flag;
2. geometric maximum-dot-product matching of transformed 600-cell vertices;
3. cell and internal-edge permutations induced from mapped vertex tuples;
4. explicit sparse Fourier isometries `Q_k`;
5. sparse products `B_k = Q_k^* (C_s Q_k)`;
6. `scipy.linalg.eigh(..., driver="evd")`, not the primary NumPy driver;
7. eigenvector residuals lifted back to the original sparse matrix.

No function implementing the primary group permutation, Fourier block or
diagonalization may be loaded.

## 3. Fixed geometry and schedules

All physical and combinatorial hypotheses remain exactly those in protocol
`2eef0e1`: `P(sd K_600)`, `19,680` internal edges, fixed boundaries,
`tau0=0.0102`, curvature-selected conserved masses, and logarithmic signed
squared-edge coordinates.

Build all schedule geometries only to obtain the common target-free stencil
catalogue.  Assemble the eleven forward representative matrices and require
their exact preregistered CSR digests.  The failed old sparse artifact and all
of its spectral fields are forbidden.

The four ambient reflection normals are obtained by SVD from the same
lexicographic base flag used in the accepted schedule-0 adversarial protocol
`b54f590`.  For column vectors use

```text
G = R_3 R_2 R_1 R_0,
G_reverse = R_0 R_1 R_2 R_3.
```

Both maps must have order 30 on vertices and cells.  On every representative
they must induce 656 cycles of length 30 on internal edges, with the reverse
edge action equal to the inverse forward action.  The maximum Euclidean
vertex-matching residual remains bounded by the already frozen `5e-8` gate.

## 4. Complete explicit Fourier census

For every representative and character `k=0,...,29`, construct the sparse
isometry whose cycle coefficients are

```text
Q_k[i_j,a] = exp(-2*pi*i*k*j/30)/sqrt(30).
```

Check `Q_k^* Q_k = I` directly.  Diagonalize in full the sixteen independent
bordered or unbordered blocks `k=0,...,15`; use the exact weights from the
primary protocol to recover dimension `19,681`.  The tangent and border
coordinate occur only at `k=0`.

For the nearest eigenpair of every block, lift the vector through `Q_k` and
measure the residual against the original unaveraged sparse matrix.  Store
every complete block spectrum, gate, residual, margin and zero-compatible
candidate.

## 5. Error rule and comparisons

The same adversarial forward-error rule frozen in `b54f590` applies without
change:

```text
|lambda| > 100 * (
    local_operator_error
  + geometric_group_average_bound
  + sparse_product_roundoff_bound
  + direct_lifted_residual
).
```

Sparse-product roundoff uses the actual maximum contribution count and

```text
gamma_q * || |Q|^T |C_s| |Q| ||_inf.
```

No primary spectrum may be parsed until all eleven adversarial spectra have
been completely constructed.  Only then compare each sorted full spectrum
with its frozen primary counterpart using 100 times the sum of both declared
forward-error bounds.

## 6. Positive and negative controls

For every representative require:

- exact source CSR digest and Lorentzian local-stencil branch;
- forward and reverse order-30 geometric edge actions with `656 x 30` cycles;
- geometric covariance inside the local operator envelope;
- all explicit `Q_k` isometries and exact weighted dimension;
- weighted trace and Frobenius/Parseval identities;
- full-spectrum agreement with the primary result after the blind phase.

Negative controls:

1. a single ambient reflection has order 2, not 30;
2. a `1e-4` corruption of the first Hessian diagonal entry must violate
   geometric covariance for every representative;
3. for the first new representative, reversing the Fourier phase without
   relabelling the sector must fail the same-sector `k=1` block comparison.

## 7. Checkpoint discipline

The only result artifact is

```text
reproducible/gravity_600cell_refined_nonhomogeneous_coxeter_census_adversarial.json
```

It may be atomically updated after each completed representative.  Normal
execution recomputes from scratch.  Explicit `--resume` may use only a prefix
checkpoint carrying the exact protocol, implementation and input hashes.
Schedule order is fixed and cannot be changed after partial results exist.

## 8. Verdict hierarchy

- Any provenance, geometry, matrix, action, covariance, isometry, invariant,
  negative-control or primary-comparison failure gives
  `ADVERSARIAL_COXETER_CENSUS_CONSTRUCTION_INVALID`; no kernel claim is
  accepted for the new representatives.
- Any zero-compatible adversarial eigenvalue gives
  `PRIMARY_COMPLETE_KERNEL_CLAIM_NOT_CORROBORATED`; report its schedule,
  character, eigenvalue, gate and lifted residual as the headline.
- If every construction control passes, every primary spectrum is reproduced
  and every adversarial eigenvalue is separated, the verdict is
  `ALL_24_SCHEDULE_KERNELS_ADVERSARIALLY_CORROBORATED` using the already
  certified time-reversal congruences.

Even the positive verdict closes only the internal-kernel question for the
fixed matrices.  It does not select a physical tick, boundary evolution,
propagation speed, `G`, or a Planck scale.  No full suite and no old sparse
census will be run.
