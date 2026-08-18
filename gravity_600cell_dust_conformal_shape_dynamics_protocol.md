# Preregistered protocol: dynamic conformal/shape closure

Date: 2026-08-18

Prior-art gate commit: `42313fb`.

No conformal/shape cross residual of `Gamma` or `Omega` has been evaluated
before this commit.  The already known global-scale invariance and the
kinematic conformal-sign result are disclosed controls, not targets.

There is one declared split per frozen matrix audit.  Both schedules, all
seven sectors and all four derivative variants are mandatory; they are not a
look-elsewhere list.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `commons/cell600.py` | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |
| `gravity_600cell_dust_conformal_supermetric.json` | `b38d55f9f575ddffd34edeaa5e835d9e10919e6d96a0c284d73c31a072675025` |
| `verify_gravity_600cell_dust_conformal_supermetric.py` | `d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4` |
| `gravity_600cell_dust_centered_jacobi.json` | `fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56` |
| `gravity_600cell_dust_centered_jacobi.npz` | `1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |

Require the inherited outcomes `CENTERED_JACOBI_CERTIFIED` and
`CONFORMAL_MAXIMAL_MINORITY_CERTIFIED`, exact hashes, `560` stored centered
arrays, sector dimensions

```text
3,2,2,2,1,1,1,
```

and the complete direct geometry controls.

## 1. Exact conformal carrier

Reconstruct each schedule's literal unsigned vertex--edge incidence matrix

```text
C[e=(u,v),w] = 1 if w=u or w=v, else 0.
```

Require the exact graph and equivariance controls inherited from the conformal
verifier, including exact row-permutation covariance between schedules.

Using the identical high-precision minimal binary-tetrahedral bases, compress
`C` in each sector.  For irrep dimension `d`, require resolved rank `r=5d` and
let `U` be its deterministic left singular basis.  With

```text
epsilon_C = 1000 eps_machine max(rows,cols) max(1,||C_d||_2),
```

require all `r` nonzero singular values above `100 epsilon_C`, all remaining
values below `10 epsilon_C`, and define the conservative conformal-projector
error

```text
eta_K = 2 epsilon_C/(s_r-2 epsilon_C)
      + 1000 eps_machine (30d).
```

The sector is open if the denominator is nonpositive.

## 2. Action-relative shape complement

For every schedule, sector and derivative variant, re-enclose the stored `M`
midpoint with its saved radii and component half-ULPs.  Put

```text
H = (M+M*)/2,
A = U* H,
S_H = ker A.
```

The complete Hermitian error is

```text
epsilon_H = ||radius_H||_F
          + 1000 eps_machine (30d) max(1,||H||_2).
```

Enclose the row operator `A` by

```text
epsilon_A = epsilon_H
          + 2 eta_K (||H||_2 + epsilon_H)
          + 1000 eps_machine (30d) max(1,||A||_2).
```

Require `A` to have resolved rank `r`.  Let `W` be its deterministic right
null basis of dimension `30d-r=25d`, and define

```text
eta_S = 2 epsilon_A/(s_r(A)-2 epsilon_A)
      + 1000 eps_machine (30d).
```

Also require the direct-sum matrix

```text
B = [U W]
```

to have all `30d` singular values resolved nonzero under the same binary SVD
floor.  Record its condition number.  No basis inside `S_H` may be rotated in
response to the dynamic matrices.

## 3. Primary invariant-subspace tests

For `X` equal to each stored `Gamma` and `Omega`, re-enclose its midpoint and
put

```text
epsilon_X = ||radius_X||_F
          + 1000 eps_machine (30d) max(1,||X||_2).
```

### Conformal invariance

The basis-free residual is

```text
R_K(X) = (I-U U*) X U.
```

Use the conservative error

```text
epsilon_K(X) = epsilon_X
             + 2 eta_K (||X||_2 + epsilon_X)
             + 1000 eps_machine (30d) max(1,||X||_2).
```

### Shape invariance

Because `S_H=ker A`, its invariant-subspace residual is

```text
R_S(X) = A X W.
```

Use

```text
epsilon_S(X)
 = epsilon_A (||X||_2 + epsilon_X)
 + ||A||_2 epsilon_X
 + (||A||_2 + epsilon_A)(||X||_2 + epsilon_X) eta_S
 + 1000 eps_machine (30d)
   max(1,||A||_2 ||X||_2).
```

For each operator norm `r_X` and corresponding error `epsilon`, classify

```text
ZERO_CONSISTENT   r_X <= 10 epsilon,
NONZERO_RESOLVED  r_X > 100 epsilon,
OPEN              otherwise.
```

The complete census contains exactly

```text
2 schedules * 7 sectors * 4 variants
* 2 operators * 2 carriers = 224
```

classifications.  Dynamic conformal/shape decoupling requires all `224` to be
`ZERO_CONSISTENT`.  One resolved nonzero classification refutes that universal
claim.

## 4. Schedule audit

For each fixed `(sector,variant,operator,carrier)` compare the scalar residual
norms between the two schedules.  Let

```text
distance = abs(r_even-r_odd),
epsilon_compare = epsilon_even + epsilon_odd
                + 1000 eps_machine (30d) max(1,r_even,r_odd).
```

Classify it as

```text
SCHEDULE_ROBUST     distance <= 10 epsilon_compare,
SCHEDULE_DEPENDENT  distance > 100 epsilon_compare,
SCHEDULE_OPEN       otherwise.
```

There are `7*4*2*2=112` comparisons.  These are diagnostics: dependence does
not turn a resolved failure in either schedule into closure.

## Frozen outcome hierarchy

1. `CONFORMAL_SHAPE_DYNAMICS_CONTROL_FAILED` for provenance, geometry,
   inherited-outcome, carrier-rank or direct-sum failure.
2. `CONFORMAL_SHAPE_CARRIER_OPEN` if either subspace error cannot be resolved.
3. `CONFORMAL_SHAPE_MIXING_REFUTED` if any required residual is
   `NONZERO_RESOLVED`.
4. `CONFORMAL_SHAPE_DYNAMICS_OPEN` if none is resolved nonzero but at least one
   is open.
5. `CONFORMAL_SHAPE_DYNAMICS_DECOUPLED` only if all `224` residuals are
   zero-consistent.

The verifier passes for every correctly resolved scientific outcome,
including refutation.

## Explicit exclusions

- no identification of `S_H` with transverse-traceless tensors or gravitons;
- no imposed Hamiltonian/diffeomorphism quotient;
- no continuum harmonic, desired degeneracy or fitted spatial operator;
- no physical time unit, dispersion, limiting speed or Planck scale;
- no refinement and no full-suite run.
