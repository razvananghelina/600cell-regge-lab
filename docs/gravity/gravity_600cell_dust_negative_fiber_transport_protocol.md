# Preregistered protocol: negative-fiber identity and tangent transport

Date: 2026-08-18

Prior-art/framing commit: `42dc0e2`.

Status: **TARGET-DISCLOSED, PREREGISTERED BEFORE ANY OLD/SHIFTED PROJECTOR OR
TANGENT-LEAKAGE NORM IS COMPUTED.**

The two rank-`15` sectors and their persistence are disclosed.  This test asks
whether their fibers, not merely their dimensions, are geometrically common
or dynamically closed.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_centered_jacobi.json` | `fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56` |
| `gravity_600cell_dust_centered_jacobi.npz` | `1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef` |
| `gravity_600cell_dust_shape_stiffness.json` | `03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868` |
| `gravity_600cell_dust_shifted_centered.json` | `265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47` |
| `gravity_600cell_dust_shifted_centered.npz` | `c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8` |
| `gravity_600cell_dust_shifted_shape_stiffness.json` | `14fe5bc91e3ae4712c6ea19b8120785e2facd364e1ceb194009123fa353a4315` |
| `gravity_600cell_dust_shifted_direct_precision.json` | `86b53f228d6cfa7326a677d881463f1b849e76bc6c9ac2b0e8aa6fd427042944` |
| `gravity_600cell_dust_two_step_full_tangent.json` | `f7fbf18535cc00dacec9a9ffa95f97f2d1847ac83073f27d39fcdb7968b0bafc` |
| `gravity_600cell_dust_two_step_full_tangent.npz` | `ce78ebf415584b1cdcf1d2cb07687135b624ad4939e0a4e54650653f7b384e6d` |
| `verify_gravity_600cell_dust_two_step_full_tangent.py` | `c1a3fb09146188c1932ab81629ab69817f2a2f19108fdf8d9e89d78b6de8f717` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_600cell_dust_conformal_supermetric.py` | `d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4` |
| `commons/cell600.py` | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |

Require the committed source outcomes, exact archive counts, the old resolved
rank and the direct shifted persistence result.  The shifted broad archive is
used deliberately: it encloses the direct midpoint and therefore gives a more
conservative projector uncertainty than the recovered sign calculation.

## Complete census

Use both parities, sectors `4,5` and all four frozen derivative variants:

```text
16 old/shifted projector pairs,
64 tangent leakage blocks = 16 pairs * {A,B,C,D}.
```

No off-target sector is inspected.

## Projector construction and error

Rebuild the exact rank-`5` conformal image and the action-selected rank-`25`
shape complement separately for the old and shifted centered `M` matrices.
Restrict `A=-V` to each shape complement and require the midpoint spectrum to
have the disclosed separated `15/10` split.

For the negative/positive gap

```text
gap = lambda_15 - lambda_14
```

and the already propagated restricted-form error `epsilon_A`, require
`gap>2 epsilon_A` and define the conservative projector uncertainty

```text
eta_eig = 2 epsilon_A/(gap-2 epsilon_A) + arithmetic_floor,
eta_P   = 2 eta_shape + eta_eig + arithmetic_floor.
```

Here `eta_shape` is the existing action-selected shape-carrier subspace bound.
Lift the negative eigenvectors back to the full `30`-position sector and form
the Hermitian projector `P`.

## Geometry-selected identity comparison

The already-derived boundary ordering supplies the identity comparison.  Set

```text
d_identity = ||P_shifted-P_old||_2,
e_identity = eta_P_shifted + eta_P_old + arithmetic_floor.
```

Per cell classify:

```text
COMMON_FIBER_RESOLVED   d_identity <= 10 e_identity,
ROTATED_FIBER_RESOLVED  d_identity > 100 e_identity,
FIBER_IDENTITY_OPEN     otherwise.
```

This is a geometric comparison, not yet dynamical transport.

## Action-selected tangent closure

Re-enclose every committed binary tangent midpoint with its stored Flint radius
plus half an ULP.  Split `T_2=[A B; C D]`.  With old/shifted projectors
`P0,P1`, compute the four fixed leakage matrices

```text
R_A = (I-P1)       A P0,
R_B = (I-P1)       B conj(P0),
R_C = (I-conj(P1)) C P0,
R_D = (I-conj(P1)) D conj(P0).
```

For each tangent block `X`, let `epsilon_X` be its complete component-radius
Frobenius bound plus the frozen arithmetic floor.  Bound the leakage error by

```text
epsilon_R = epsilon_X
          + (eta_P0 + eta_P1 + eta_P0 eta_P1)
            (||X||_2 + epsilon_X)
          + arithmetic_floor.
```

Classify each block:

```text
LEAKAGE_ZERO_CONSISTENT  ||R||_2 <= 10 epsilon_R,
LEAKAGE_NONZERO_RESOLVED ||R||_2 > 100 epsilon_R,
LEAKAGE_OPEN             otherwise.
```

## Frozen outcome hierarchy

Use the first applicable branch:

1. `NEGATIVE_FIBER_TRANSPORT_CONTROL_FAILED` if provenance, carrier ranks,
   gaps, projector errors, tangent finiteness or census completeness fail;
2. `NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED` if any of the `64` leakage blocks
   is nonzero-resolved;
3. `NEGATIVE_FIBER_TANGENT_CLOSURE_OPEN` if none is nonzero-resolved but any is
   open;
4. `NEGATIVE_FIBER_TANGENT_CLOSURE_CERTIFIED` only if all `64` blocks are
   zero-consistent.

Report the independent identity-fiber classification counts regardless of the
tangent outcome.  A tangent refutation is not repaired by a small identity
distance.

## Scope and exclusions

- no Procrustes, polar or permutation optimization;
- no fitted phase lift or omission of an inconvenient `A,B,C,D` block;
- no claim that failure excludes every possible constraint-derived lift;
- no physical instability, graviton, wave, mass, speed or continuum claim;
- no full-suite run.
