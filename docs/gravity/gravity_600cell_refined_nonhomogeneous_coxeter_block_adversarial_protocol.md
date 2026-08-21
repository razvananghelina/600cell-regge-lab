# Adversarial protocol: geometric Coxeter replication of schedule 0

Date: 2026-08-21

Status: frozen after the primary result `fd75a55` and before constructing any
geometric reflection or adversarial Fourier block.

## 1. Claim under attack

The primary verifier reports, for the fixed lexicographic refined slab and
its bordered internal Hessian,

```text
ker(C0) = span(n0)
```

with a least zero-exclusion margin ratio near 6.  This remains **PRIMARY
ONLY** until the present mechanically different route corroborates it.

All physical hypotheses are unchanged: `P(sd K_600)`, `tau0=0.0102`, the
curvature-selected conserved vertex masses, logarithmic signed-squared-edge
coordinates, fixed spatial boundaries, and all `19,680` internal directions.

## 2. Independent geometric group action

The primary construction obtained the Coxeter action as the coloured-graph
automorphism commuting with the four right chamber adjacencies.  The
adversarial construction may not use that permutation.

Instead, assign to each coarse cell the normalized Euclidean sum of its
600-cell vertex vectors in `R4`.  The four centres in the lexicographic base
flag form a spherical Coxeter chamber.  For rank `r`, compute the unit normal
to the span of the other three centres by an SVD and construct the ambient
reflection

```text
R_r = I - 2 n_r n_r^T.
```

For column vectors, freeze the sequential word `(0,1,2,3)` as

```text
G = R_3 R_2 R_1 R_0.
```

Map every transformed coarse vertex to the unique stored 600-cell vertex of
maximum dot product.  The map must be bijective and its maximum Euclidean
matching residual must not exceed `5e-8`, a bound fixed from the `1e-10`
coordinate rounding in `build_600cell`, before observing the result.

Map every higher-rank coarse cell by mapping its vertex tuple, then induce
the two-layer internal-edge permutation.  It must preserve rank and the
schedule, have order 30, and consist of 656 cycles of length 30.  No chamber
BFS or colour-centralizer permutation from the primary verifier is loaded.

As a convention attack, the reverse matrix word

```text
G_reverse = R_0 R_1 R_2 R_3
```

must also induce an order-30 carrier symmetry.  Its complete spectrum must
agree after sector relabelling; at minimum the explicitly constructed
nontrivial sector pair is compared.

## 3. Independent block construction

The primary verifier accumulated each block entry directly from COO matrix
entries and phase differences.  That decisive step is forbidden here.

For every cyclic character construct the explicit sparse isometry

```text
Q_k in C^(19680 x 656),
Q_k[i_j,a] = exp(-2*pi*i*k*j/30)/sqrt(30)
```

with one nonzero per row.  Verify `Q_k^* Q_k = I` directly.  Form each block
only through sparse matrix products

```text
B_k = Q_k^* (C0 Q_k).
```

The bordered invariant block uses `Q_0^* n0`.  All independent sectors
`k=0,...,15` are diagonalized completely with `scipy.linalg.eigh` using the
divide-and-conquer driver `evd`, rather than the primary NumPy call.  Every
eigenpair residual is checked against the original sparse matrix after
lifting the nearest vector back with `Q_k`; block residuals alone are not the
only evidence.

## 4. Error and comparison gates

The primary source matrix must reproduce its frozen CSR digest.  Its measured
geometric covariance defect is recomputed using the new permutation and
added to the same local operator error.  Sparse Fourier multiplication
roundoff is bounded conservatively by its measured Hermitian defect plus

```text
gamma_q * || |Q|^T |C0| |Q| ||_inf,
```

where `q` is the maximum actual number of scalar contributions to a block
entry and `gamma_q=q*u/(1-q*u)`.

An adversarial eigenvalue is separated from zero only when

```text
|lambda| > 100 * (
    local_operator_error
  + geometric_group_average_bound
  + sparse_product_roundoff_bound
  + direct_lifted_residual
).
```

After all blocks are complete, reconstruct both full weighted spectra and
compare the sorted eigenvalue lists from the primary and adversarial
artifacts.  The maximum difference must be below 100 times the sum of the two
declared forward-error bounds.  The primary list is not read until the
adversarial spectrum has been constructed.

## 5. Controls and verdict hierarchy

Positive controls:

1. the geometric reflections reproduce the `H4=(3,3,5)` pair orders from
   their matrices and permute all coarse vertices/cells bijectively;
2. the internal action has `656 x 30` cycles;
3. every explicit `Q_k` is an isometry and the weighted dimensions sum to
   `19,681` after bordering;
4. weighted trace and Frobenius/Parseval invariants reproduce the full
   bordered matrix;
5. the full adversarial spectrum agrees with the frozen primary spectrum.

Negative controls:

1. replacing `G` by a single chamber reflection must give order 2, not 30;
2. adding `1e-4` to the first diagonal Hessian entry must violate geometric
   Coxeter covariance by more than the frozen gate;
3. reversing one Fourier phase sign without relabelling `k` must exchange a
   non-real sector with its conjugate and fail the same-sector block test.

Verdicts:

- if all construction controls pass, both full spectra agree and every
  adversarial eigenvalue is separated, the schedule-0 claim is
  `ADVERSARIALLY_CORROBORATED`;
- if a zero-compatible sector appears, the primary positive result is
  refuted or unresolved, and that sector is the headline;
- if the geometric carrier, isometries, controls or full-spectrum comparison
  fail, the outcome is `ADVERSARIAL_CONSTRUCTION_INVALID`; no kernel claim is
  accepted.

No conclusion about the other schedules, propagation, gravitons, a selected
tick, `c`, `G` or Planck units is permitted.  No full suite is run.
