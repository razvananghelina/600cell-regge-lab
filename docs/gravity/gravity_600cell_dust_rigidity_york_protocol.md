# Preregistered protocol: dynamic closure of the 600-cell self-stress carrier

Date: 2026-08-18

Prior-art gate commit: `a318d6e`.

No rigidity-sector decomposition of a centered matrix, cross-block norm,
self-stress leakage or self-stress restricted spectrum has been inspected
before this commit.  The theorem-predicted global ranks `470` and `250` are
disclosed controls, not targets discovered by the calculation.

There is one canonical embedded rigidity map and one self-stress carrier.
Schedules, sectors and derivative variants are audits, not alternative
candidates.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `commons/cell600.py` | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |
| `gravity_600cell_dust_conformal_supermetric.json` | `b38d55f9f575ddffd34edeaa5e835d9e10919e6d96a0c284d73c31a072675025` |
| `verify_gravity_600cell_dust_conformal_supermetric.py` | `d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4` |
| `gravity_600cell_dust_centered_jacobi.json` | `fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56` |
| `gravity_600cell_dust_centered_jacobi.npz` | `1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |

Require the centered and conformal outcomes respectively
`CENTERED_JACOBI_CERTIFIED` and
`CONFORMAL_MAXIMAL_MINORITY_CERTIFIED`, exact recorded hashes, 560 centered
arrays, and the two schedules with dimensions `3,2,2,2,1,1,1`.

## 1. Literal rigidity controls before dynamics

Reconstruct the unit 600-cell vertices `x_v` and each schedule's literal
edge order.  With the common squared edge length `L2`, form

```text
R[e=(u,v),(u,a)] =  2 (x_u-x_v)_a / L2,
R[e=(u,v),(v,a)] = -2 (x_u-x_v)_a / L2.
```

All other entries vanish.  Form the radial injection

```text
J[(v,a),v] = (x_v)_a
```

and the block tangent projector

```text
P_T[v] = I_4 - x_v x_v^T.
```

Require, within the complete fixed-coordinate arithmetic envelope:

```text
R J = C,
rank R = 470,
rank(R P_T) = 354,
rank(C, R P_T) = 470,
dim(im C intersect im(R P_T)) = 4,
dim ker R^T = 250.
```

The first identity fixes the normalization.  The ranks are controls from
Whiteley's theorem, not scientific outcomes.  Also require exact row
permutation covariance between the two schedule edge orders.

For a matrix `A` use

```text
epsilon_A = 1000 eps_machine max(rows,cols) max(1,||A||_2).
```

Singular values above `100 epsilon_A` are nonzero-resolved, below
`10 epsilon_A` are zero-consistent and the rest open.  Any failure of the
theorem-predicted ranks is a carrier-control failure, not evidence against
the theorem.

## 2. Minimal-sector self-stress carriers

Use the identical high-precision minimal bases and edge order as the frozen
centered matrices.  For each irrep dimension `d` form

```text
R_d = (I_30 tensor B_d)* R.
```

Do not preregister individual sector ranks.  Resolve them blindly and require
only the theorem-controlled restoration identities

```text
sum_d d rank(R_d) = 470,
sum_d d (30d-rank(R_d)) = 250.
```

Let `U_d` span `im R_d` and `W_d` span its Euclidean edge-orthogonal
complement.  Thus `W_d` is the minimal-sector self-stress carrier.  Require
projector residuals below `10 epsilon_R`.

For the last resolved rigidity singular value `s_gap`, define the conservative
subspace envelope

```text
eta_S = 2 epsilon_R/(s_gap-2 epsilon_R)
      + 1000 eps_machine (30d)
```

when `s_gap>2 epsilon_R`; otherwise the sector is open.

## 3. Primary dynamic-decoupling test

For every schedule, sector, variant and `X` in `M,N,V`, re-enclose the stored
matrix midpoint with its radius plus component half-ULPs.  Define

```text
L_X = U_d* X W_d,
R_X = W_d* X U_d.
```

Both must vanish for a block-decoupled bilinear recurrence.  Use

```text
epsilon_X = ||radius_X||_F
          + 1000 eps_machine (30d) max(1,||X||_2),
epsilon_cross = epsilon_X + 2 eta_S ||X||_2.
```

Classify each operator norm as

```text
ZERO_CONSISTENT   norm <= 10 epsilon_cross,
NONZERO_RESOLVED  norm > 100 epsilon_cross,
OPEN              otherwise.
```

No post-result rotation of `U_d` or `W_d` is allowed.

For `Gamma=M^-1N` and `Omega=M^-1V`, record the normalized recurrence
leakages

```text
||(I-W_d W_d*) Gamma W_d||_2,
||(I-W_d W_d*) Omega W_d||_2,
```

with the same classification and their stored matrix envelopes.  These are
equivalent to the `U_d` cross blocks because `U_d,W_d` exhaust the minimal
edge sector.

The self-stress carrier is dynamically closed only if every left and right
cross block for `M,N,V` and both normalized leakages is zero-consistent in
all `2*7*4` audits.

## 4. Schedule comparison

For `Gamma` and `Omega`, compare the ordered singular spectra of the
restricted blocks

```text
W_d* X W_d
```

between schedules.  Use the sum of both matrix/subspace envelopes plus the
binary SVD floor and classify as `SCHEDULE_ROBUST`, `SCHEDULE_DEPENDENT` or
`SCHEDULE_OPEN` with the usual `10/100` bands.  These comparisons are
reported diagnostics and do not override the closure result.

## Frozen outcome hierarchy

1. `RIGIDITY_YORK_CONTROL_FAILED` for provenance, geometry, normalization,
   row-permutation, high-precision basis or theorem-rank failure.
2. `RIGIDITY_YORK_CARRIER_OPEN` if a minimal carrier rank/subspace is open.
3. `RIGIDITY_YORK_DECOUPLING_REFUTED` if any required cross block or
   normalized leakage is resolved nonzero.
4. `RIGIDITY_YORK_DECOUPLING_OPEN` if none is resolved nonzero but at least
   one is open.
5. `RIGIDITY_YORK_DECOUPLING_CERTIFIED` only if every required quantity is
   zero-consistent.

The verifier passes for any correctly resolved scientific outcome, including
refutation.

## Explicit exclusions

- `250` is not called a new count, a physical quotient or a graviton count;
- no tangent-vertex direction is discarded as gauge by assumption;
- no desired two-polarization spectrum, continuum harmonic or wave speed;
- no fitted basis, coefficient or threshold;
- no proper time, Planck scale, matter mass or Standard-Model target;
- no refinement and no full-suite run.
