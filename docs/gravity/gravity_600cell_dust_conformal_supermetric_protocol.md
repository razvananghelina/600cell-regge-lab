# Preregistered protocol: canonical conformal restriction of the Regge kinetic form

Date: 2026-08-18

Prior-art gate commit: `96dd1ff`.

No stored centered-matrix eigenvector, conformal restriction, principal
angle or invariance residual has been inspected before this commit.  The
previously committed target-free fact that the full Hermitian inertia is
`120:600` is disclosed.

There is exactly one candidate map.  The two schedules, seven sectors and
four derivative variants are repeated audits of that map, not 56 candidates.
No normalization of the map can change any outcome below.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_centered_jacobi.json` | `fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56` |
| `gravity_600cell_dust_centered_jacobi.npz` | `1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef` |
| `verify_gravity_600cell_dust_centered_jacobi.py` | `359b8d7642746c2dc22e304353e3b83104874badd86755de4f8f9e6f25e56a20` |
| `gravity_600cell_dust_full_boundary_tangent.json` | `4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |

Require the centered outcome `CENTERED_JACOBI_CERTIFIED`, `560` source
arrays, the recorded source archive hash, both schedules and sector
dimensions `3,2,2,2,1,1,1`.

## 1. Literal canonical map

For each schedule reconstruct the old boundary edge order directly from the
fixed geometry and the same lexicographically ordered binary-tetrahedral
actions used to construct the source blocks.  Form the integer matrix

```text
C[e=(u,v), w] = 1 when w=u or w=v, and 0 otherwise.
```

Require exactly:

- shape `720 x 120`;
- two ones in every edge row and twelve ones in every vertex column;
- `C^T C = 12 I + A`, where `A` is the literal 600-cell adjacency matrix;
- a connected adjacency graph containing a triangle;
- exact injectivity by the odd-cycle argument, hence rank `120`;
- exact equivariance for all 24 group elements;
- the same literal old-boundary edge set for both schedules, with each
  schedule's independently reconstructed orbit ordering a bijection of that
  set;
- an exact row permutation relating the two incidence matrices.  Each
  centered matrix must be paired with the incidence matrix in its own
  schedule-specific orbit order; equality of the two row sequences is not
  required.

These are controls, not data-dependent outcomes.

Protocol correction: the first committed version incorrectly required the
two orbit-ordered row sequences to be literally identical.  A geometry-only
check performed before any centered eigenspace was read found the first
difference at row `72`, although both sequences are bijections of the same
720 literal edges and are related by an exact permutation.  The scientific
test is invariant under that paired permutation.  This correction is being
committed separately rather than hidden by rewriting the original
preregistration.

## 2. Minimal-sector conformal images

Reconstruct the same seven 100-decimal minimal regular-representation bases
`B_d` as the full-boundary tangent verifier.  In the old-edge order define

```text
E_d = I_30 tensor B_d,
C_d = E_d* C.
```

The full vertex carrier consists of five free regular group orbits, so the
preregistered sector rank is `5d`.  Calculate an orthonormal column basis
`U_d` for `im C_d` by SVD.

Define the geometry arithmetic envelope

```text
epsilon_C = 100 eps_machine max(rows,cols) max(1, ||C_d||_2).
```

A singular value is nonzero-resolved above `100 epsilon_C`, zero-consistent
below `10 epsilon_C`, and open otherwise.  Require exactly `5d` resolved
nonzero singular values and all remaining singular values zero-consistent.
Also require `||(I-U_d U_d*) C_d||_2 <= 10 epsilon_C` and the high-precision
group-basis controls below `1e-70`.

Sector ordering is matched to the frozen source by the complete tuple

```text
(irrep dimension, old central eigenvalue, splitter group index),
```

not by dimension alone.

## 3. Primary restricted-inertia gate

For every schedule, sector and all four derivative variants, re-enclose the
stored `M` midpoint using its stored radius plus half an ULP per component,
then form

```text
H = (M+M*)/2,
G = U_d* H U_d.
```

Use a conservative operator envelope

```text
epsilon_H = ||R_H||_F
          + 1000 eps_machine n max(1,||H_mid||_2),
epsilon_G = epsilon_H
          + 1000 eps_machine n max(1,||H_mid||_2).
```

Here `R_H[i,j]=(R_M[i,j]+R_M[j,i])/2`, including the half-ULP
re-enclosure, and `n=30d`.  The second term covers binary projection and
Hermitian eigensolve arithmetic.  Classify eigenvalues of both `H` and `G`
using the already frozen bands:

```text
POSITIVE_RESOLVED  lambda >  100 epsilon,
NEGATIVE_RESOLVED  lambda < -100 epsilon,
ZERO_CONSISTENT   |lambda| <  10 epsilon,
OPEN               otherwise.
```

Every full `H` must independently reproduce `(5d,25d,0,0)` and every
restriction `G` must have `(5d,0,0,0)`.  If so, the literal conformal image
is a maximal positive subspace in every audit.  Since positive is the
committed minority sign in the stored action convention, this is invariantly
the maximal-minority statement up to reversal of the whole action.

Record all restricted eigenvalue extrema, margins in error units,
conditions and schedule/variant distances.  Do not fit or rotate `C`.

## 4. Secondary structural diagnostics

Only after the primary restriction is computed, use the operational-primary
midpoint to form the Euclidean positive spectral basis `W_d` of `H` and
report

```text
projector_distance = ||U_d U_d* - W_d W_d*||_2,
spectral_angle_max  = asin(min(1,projector_distance)),
leakage             = ||(I-U_d U_d*) H U_d||_2.
```

Let `gap` be the distance between the nearest positive and negative
eigenvalues of `H`.  If `gap <= 2 epsilon_H`, the diagnostic is open.
Otherwise use

```text
epsilon_P = 2 epsilon_H / (gap - 2 epsilon_H)
          + 1000 eps_machine n,
epsilon_L = epsilon_H
          + 1000 eps_machine n max(1,||H||_2).
```

Classify projector equality and invariance independently:

```text
IDENTIFIED  value <= 10 epsilon,
SEPARATED   value > 100 epsilon,
OPEN        otherwise.
```

These labels are **STRUCTURAL** and never alter the primary outcome.  A
maximal definite subspace need not equal a spectral subspace defined using
an auxiliary Euclidean norm.

## Frozen outcome hierarchy

1. `CONFORMAL_SUPERMETRIC_CONTROL_FAILED` for any provenance, geometry,
   carrier, equivariance, sector-match or source-cardinality failure.
2. `CONFORMAL_SUPERMETRIC_SECTOR_OPEN` if a sector rank is not resolved.
3. `CONFORMAL_SUPERMETRIC_RESTRICTION_OPEN` if any full or restricted sign
   is open and none has already refuted the claim.
4. `CONFORMAL_MAXIMAL_MINORITY_REFUTED` if a restricted matrix has a resolved
   negative direction, or a full matrix does not reproduce the required
   inertia.
5. `CONFORMAL_MAXIMAL_MINORITY_CERTIFIED` only if all 56 full matrices have
   `(5d,25d,0,0)` and all 56 restrictions have `(5d,0,0,0)`.

The spectral-subspace and invariance labels are reported separately.

## Explicit exclusions

- no desired continuum tensor spectrum or polarization count;
- no Whitney/Hodge, graph, Kähler--Dirac or Lichnerowicz target;
- no constraint quotient or graviton label;
- no fitted Schur coefficient or post-result rotation;
- no proper-time normalization, dispersion, `c`, Planck scale or particle
  mass;
- no refinement and no full-suite run.
