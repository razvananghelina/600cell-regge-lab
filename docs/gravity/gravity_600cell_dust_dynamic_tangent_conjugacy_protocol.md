# Preregistration: finite geometric conjugacy census

Date: 2026-08-17

Prior-art gate: `0129053`.

Status: frozen before enumerating a cross-parity carrier map or evaluating an
intertwining residual.

## 1. Frozen sources

Require these SHA-256 values:

```text
one-slab carrier/action source
ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf

two-slab gluing source
9ea55dab1fd2f4e9ee643247f5d35599c5894cf77970fc2006fe3d8ac22edf37

two-slab gluing artifact
a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77

blind dynamic tangent artifact
1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5.
```

The carrier source must retain `43/43` checks.  The two models must have the
same 120 spatial vertices, 720 spatial edges, five cover cells, order-24
ordered stabilizer and `30+35+30` orbit dimensions.  Their schedules must be
the even ordering and the ordering obtained by swapping only phases zero and
one.

## 2. Stage A: enumerate before comparing

Stage A may load the carrier and gluing data but must not load either tangent
matrix, eigenvalue, trace or spectral target.  It writes and commits the full
candidate list before Stage B exists or runs.

### 2.1 Complete-slab maps

For every one of the finite `H4` vertex actions, test exactly two extensions to
the two-layer vertex set:

```text
direct:   (v,t) -> (g(v),t),
reversed: (v,t) -> (g(v),1-t).
```

Retain an action only if it maps the complete set of 2400 even four-simplices
onto the complete set of 2400 odd four-simplices.  Derive, rather than choose,
its two endpoint orbit permutations:

- direct: `Q_old : even old -> odd old` and
  `Q_final : even final -> odd final`;
- reversed: `Q_old_final : even old -> odd final` and
  `Q_final_old : even final -> odd old`.

Deduplicate identical permutation pairs while recording their action
multiplicity.  These are the **canonical geometric** candidates.

### 2.2 Boundary-only maps

Independently enumerate every distinct permutation of the 30 old-boundary
orbits induced by an `H4` action that maps the entire even old-orbit partition
onto the odd partition.  Also include the unique map obtained by identifying
the same physical edge sets in the two partitions.  Deduplicate and record all
provenance labels and action multiplicities.

These maps are **STRUCTURAL**, not complete-slab symmetries.  Report before
comparison:

```text
N_direct_slab,
N_reversed_slab,
N_boundary,
the complete ordered multiset of permutation records.
```

The enumeration artifact must explicitly record
`tangent_matrices_parsed=false` and `spectral_target_parsed=false` and be
committed before Stage B.

## 3. Stage B: frozen canonical lifts

Only after the enumeration artifact is committed may Stage B load the two
stored `60 x 60` tangent matrices.  If a 30-permutation `q` maps source indices
to target indices, define `Q[target,source]=1` and only the two lifts

```text
C(Q)=diag(Q,Q),       C^T Omega C = Omega,
K(Q)=diag(Q,-Q),      K^T Omega K = -Omega.
```

No dense, diagonal, Schur-block or eigenvector-derived intertwiner is allowed.

For every complete direct candidate test

```text
T_odd C(Q_old) = C(Q_final) T_even.
```

For every complete reversed candidate test the covariant time-reversal
identity

```text
T_odd K(Q_final_old) T_even = K(Q_old_final).
```

For every boundary-only `Q`, test both deliberately weaker identities

```text
direct:   T_odd C(Q) = C(Q) T_even,
reversed: T_odd K(Q) T_even = K(Q).
```

Thus the boundary look-elsewhere denominator is exactly `2*N_boundary`.

## 4. Frozen calibration

Read each stored tangent-map calibration as

```text
delta_q = epsilon_t_q + 60*0.5e-50.
```

Use Frobenius residuals.  For a direct identity set

```text
u_direct = sqrt(60)*(delta_even+delta_odd).
```

For a reversed product, let `sigma_q` be the stored largest full-map singular
value plus its `epsilon_svd` and set

```text
u_reversed = sqrt(60)*(
    delta_odd*sigma_even
  + sigma_odd*delta_even
  + delta_odd*delta_even).
```

Classify each frozen candidate:

- `PASS` when residual `<=10*u`;
- `FAIL` when residual `>100*u`;
- `OPEN` otherwise.

Report every residual, uncertainty, ratio and the exact PASS/total fractions.
No candidate or relation may be added after inspection.

## 5. Mechanical verdicts

Complete-slab verdict:

- any PASS: `COMPLETE_SLAB_COVARIANCE_DERIVED`;
- candidates exist, none PASS, none OPEN:
  `COMPLETE_SLAB_COVARIANCE_REFUTED`;
- no candidates: `NO_CROSS_PARITY_H4_SLAB_ISOMORPHISM`;
- otherwise: `COMPLETE_SLAB_COVARIANCE_OPEN`.

Boundary-only verdict:

- any PASS: `BOUNDARY_INTERTWINER_STRUCTURAL`;
- none PASS, none OPEN: `BOUNDARY_INTERTWINER_REFUTED`;
- otherwise: `BOUNDARY_INTERTWINER_OPEN`.

Combined outcome:

- complete-slab PASS:
  `DYNAMIC_TANGENT_GEOMETRIC_COVARIANCE_DERIVED`;
- no complete-slab PASS but boundary PASS:
  `DYNAMIC_TANGENT_BOUNDARY_COVARIANCE_ONLY`;
- any unresolved candidate and no PASS:
  `DYNAMIC_TANGENT_COVARIANCE_OPEN`;
- otherwise:
  `DYNAMIC_TANGENT_ISOSPECTRALITY_UNEXPLAINED`.

Even a complete reversed covariance is not automatically a similarity unless
its two endpoint permutations coincide.  Exact isospectrality remains a
separate theorem.

## 6. Scope and execution

The result concerns only the finite order-24 quotient and the first accepted
dynamic dust tick.  It does not establish nonlinear schedule independence,
refinement, full 720-edge dynamics, gravitons, a causal cone or a physical
clock.

Register each verifier before its first execution.  Run only the two targeted
verifiers; do not run the full suite.
