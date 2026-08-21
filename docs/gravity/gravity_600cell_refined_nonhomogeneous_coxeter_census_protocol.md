# Protocol: exhaustive Coxeter census for all schedule representatives

Date: 2026-08-21

Status: frozen after the accepted schedule-0 result `66e47a7` and before
assembling or diagonalizing any Coxeter block for schedules other than 0.

## 1. Prior art and purpose

The prior-art gate `295e90d` remains controlling.  Finite-group block
diagonalization is standard representation theory; it is not new physics.
The only instance-specific question is whether the complete internal Hessian
has any null direction beyond the analytic product-duration line.

Schedule 0 has already been exhaustively computed twice by mechanically
different routes.  It is retained here only as a frozen positive and drift
control.  The new target-free question is

```text
ker(C_s) = span(n_s)
```

for each of the other eleven representatives of the twelve time-reversal
pairs.  No continuum, propagation, particle, tick, `c`, `G`, or Planck-scale
target may be loaded or compared.

## 2. Fixed hypotheses and schedule census

Every matrix uses exactly:

```text
carrier                    P(sd K_600)
spatial f-vector           (2640, 17040, 28800, 14400)
slab edges                 53760
boundary edges             34080
internal edges             19680
input duration             tau0 = 0.0102
vertex masses              m_v = K_v/(8*pi)
coordinates                logarithmic signed squared edge lengths
boundary condition          both spatial boundaries fixed
```

The 24 labelled schedules are the lexicographic permutations of `(0,1,2,3)`.
Time reversal maps an order to its reversed tuple.  The twelve representative
indices and their reverse partners are frozen as

```text
(0,23), (1,17), (2,21), (3,11), (4,15), (5,9),
(6,22), (7,16), (8,19), (10,13), (12,20), (14,18).
```

The machine-readable structural manifest is

```text
reproducible/gravity_600cell_refined_nonhomogeneous_csr_manifest.json
```

It contains only the schedule orders, pair map and frozen CSR digests.  It
was extracted from the failed old sparse artifact whose SHA-256 is
`4a05968c68f8e6a35a1308ddf6114bb19b7106f214bfdcf798e7af2387bddec1`.
That old artifact's `LOCAL_EXTENSION_INVALID` verdict is binding: none of its
eigenvalue or kernel fields may be used as evidence, as targets, or as solver
initialization.  The manifest deliberately contains no spectral field.

## 3. Source-matrix reconstruction

Build all 24 local slab geometries and one shared high-precision catalogue of
all local angle and area stencils.  Assemble every complete `19680 x 19680`
real-symmetric internal Hessian `C_s`; do not assemble only an orbit quotient.
Each matrix must reproduce its preregistered CSR digest in the structural
manifest.

For every schedule, independently check:

1. stationarity of all local internal gradients;
2. reality and reciprocity of the raw Hessian;
3. the analytic relation `C_s n_s = 0` inside the local operator and sparse
   multiplication envelope;
4. congruence with the explicitly assembled reverse schedule under the
   mechanically derived layer-reversal permutation.

Only the twelve forward representatives are diagonalized.  A reverse member
inherits its kernel statement only after the full-matrix congruence check for
that pair passes.

## 4. Fixed Coxeter action and exhaustive blocks

Use the colour-preserving **left** Coxeter action from the primary route,
obtained by propagating the image of the base chamber under the frozen word
`(0,1,2,3)`.  Construct it once, but induce and verify its internal-edge
permutation separately for every representative schedule.

Before any diagonalization, each permutation must:

- preserve every internal edge of its schedule;
- have order exactly 30;
- consist of exactly 656 cycles of length 30;
- commute with `C_s` inside the measured group-averaging envelope.

For every representative and every cyclic character `k=0,...,29`, use the
same normalized Fourier convention already frozen in `acfe795`.  Assemble
all sixteen independent Hermitian blocks `k=0,...,15`; sectors `1,...,14`
have spectral weight two and sectors 0 and 15 weight one.  The tangent and
border coordinate occur only in the invariant sector, producing the exact
weighted dimension

```text
657 + 14*(2*656) + 656 = 19681.
```

Diagonalize every independent block completely.  A truncated Ritz window is
forbidden.  Verify weighted trace and Frobenius/Parseval identities and build
sector 29 explicitly to check it against the conjugate of sector 1.

## 5. Frozen numerical error rule

For all powers `r=0,...,29` measure

```text
delta_r = ||C_s - P_r^T C_s P_r||_inf
```

and use their mean as the group-averaging bound.  Bound Fourier assembly by
the actual contribution counts and absolute entry sums.  For every computed
eigenpair use its direct Hermitian residual.  An eigenvalue is separated from
zero only if

```text
|lambda| > 100 * (
    local_operator_error
  + group_averaging_bound
  + Fourier_assembly_roundoff_bound
  + direct_eigenpair_residual
).
```

No tolerance may be changed after seeing any spectrum.  Store every complete
block spectrum, gate, residual, minimum margin and zero-compatible candidate.

## 6. Preregistered controls

1. Schedule 0 must reproduce the frozen complete primary spectrum with a
   maximum discrepancy below the sum of the two declared forward-error
   envelopes.  The frozen artifact is read only after the new schedule-0
   spectrum has been constructed.
2. The forbidden right-product chamber convention must fail to descend to
   rank cells or fail covariance, as in the frozen schedule-0 control.
3. Adding `1e-4` to the first Hessian diagonal entry of each representative
   must exceed that representative's covariance gate.
4. Omitting any cyclic sector or assigning the wrong spectral weight must
   fail the exact dimension, trace or Parseval census.
5. Every reverse pair must pass full-matrix congruence; equality of only a
   reduced orbit matrix is insufficient.

## 7. Verdict hierarchy

- If any source digest, topology, reversal, action, covariance, Fourier,
  dimension, trace, Parseval or control check fails, the result is
  `COXETER_CENSUS_CONSTRUCTION_INVALID`; no kernel statement is accepted.
- If construction passes but any eigenvalue is compatible with zero, report
  its schedule, order, character, eigenvalue, gate and lifted residual.  The
  result is `ADDITIONAL_INTERNAL_KERNEL_CANDIDATE_FOUND`; this is the
  headline, even if the other schedules pass.
- If every eigenvalue of every bordered representative is separated, the
  primary result is `ALL_12_REPRESENTATIVE_KERNELS_ARE_PRODUCT_DURATION_LINES_PRIMARY`.
  Schedule 0 remains adversarially accepted; the other eleven remain
  **PRIMARY ONLY** pending a mechanically different geometric-reflection and
  sparse-`Q` replication.

Even the final positive verdict would not select a physical tick.  It would
only say that, at the fixed input `tau0`, the product-duration tangent is the
unique internal degeneracy for every labelled schedule up to time reversal.
The duration-selection and boundary-evolution question is separate.

## 8. Execution and artifacts

The registered verifier will be

```text
reproducible/verify_gravity_600cell_refined_nonhomogeneous_coxeter_census.py
```

and its only result artifact will be

```text
reproducible/gravity_600cell_refined_nonhomogeneous_coxeter_census.json
```

No full suite and no old sparse factorization census will be run.  The
implementation must be committed and registered before its first execution.
