# Protocol: action-weighted longitudinal identity of the negative shape carrier

Date: 2026-08-19

Prior-art/framing commit: `9bd990a`.

This is a target-disclosed test.  The already known rank-`15` negative spaces
in sectors `4,5` are compared with a geometry/action carrier fixed without a
spectral fit.  No continuum harmonic, polarization, speed or mass target is
loaded.

## 1. Frozen inputs

Require byte-exact provenance for:

```text
docs/gravity/gravity_600cell_dust_action_york_negative_prior_art.md
SHA-256 e5728865b8498c5750cdbf45d9d93938c530ba40d95c441452f34c02bf00cd1d

reproducible/gravity_600cell_dust_centered_jacobi.json
SHA-256 fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56

reproducible/gravity_600cell_dust_centered_jacobi.npz
SHA-256 1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef

reproducible/gravity_600cell_dust_shape_stiffness.json
SHA-256 03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868

reproducible/gravity_600cell_dust_rigidity_york.json
SHA-256 251851c08f81ba2f0d41c2d0da428ab11f1ba918b9cb59e0a1e347143c883981

reproducible/verify_gravity_600cell_dust_shape_stiffness.py
SHA-256 d4f0a9a805910de37011ba70f407907daa2d11c650aeea22e571ab867282a44c

reproducible/verify_gravity_600cell_dust_rigidity_york.py
SHA-256 deba8d9f9bca4a5848134943ec77544e5487d44a59c44234f632b6f2aeb51382

reproducible/verify_gravity_600cell_dust_negative_shape_dynamics.py
SHA-256 6e7659ca398037e806f9a35a9f3db3d6035f992a8655699b47a2519b0c37453e

commons/cell600.py
SHA-256 ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

The upstream outcomes must remain:

```text
CENTERED_JACOBI_CERTIFIED,
SHAPE_STIFFNESS_NEGATIVE_MODES_RESOLVED,
RIGIDITY_YORK_DECOUPLING_REFUTED.
```

## 2. Complete geometry and symmetry reconstruction

Reconstruct the literal 600-cell, both schedule edge orders, the conformal
incidence `C`, the full rigidity differential `R`, the tangent projector

```text
P_tan,v=I_4-x_v x_v*,       D=R P_tan,
```

and all seven frozen minimal binary-tetrahedral sectors.  Require the known
global controls

```text
rank C=120,
rank R=470,
rank D=354,
dim(im C intersection im D)=4.
```

No selected stiffness matrix may be read before these carrier controls pass.

## 3. Action-weighted conformal/shape split

For every

```text
2 schedules * 7 sectors * 4 derivative variants = 56 cells,
```

reconstruct the midpoint/radius pairs `M,V` from the centered archive.  Let

```text
H=(M+M*)/2,
K=im C,
S_H=ker(C*H).
```

Use Euclidean-orthonormal bases `U,W` of `K,S_H`.  Require the inherited
rank, direct-sum and definite-kinetic gates.  The oblique projector along `K`
onto `S_H` is reconstructed as

```text
Q=[U,W],
P_S=W [0,I] Q^-1.
```

Check directly that

```text
P_S^2=P_S,       P_S U=0,       U*H P_S=0
```

inside the complete propagated error.

## 4. Longitudinal carrier before the target comparison

For all seven sectors form

```text
L_H=im(P_S D),
```

expressed in `W` coordinates.  Record, before selecting sectors `4,5` in the
output logic:

- `dim L_H`;
- `dim T_H=dim S_H-dim L_H`, where
  `T_H=ker(L_H* B)` and `B=-W*H W>0`;
- the rank gap and complete Wedin-style subspace error;
- the global multiplicity-restored sums.

Require

```text
sum d*dim(L_H,d)=350,
sum d*dim(T_H,d)=250.
```

The complete error includes source matrix radii, the inherited conformal and
shape subspace errors, the conditioning of `Q`, the fixed 600-cell coordinate
closure envelope and a `1000 eps_machine` arithmetic term.  Labels use the
existing bands: at most `10*error` is consistent with zero/equality; more than
`100*error` is resolved nonzero/separated; the interval between is `OPEN`.

## 5. Target-disclosed comparison

Only after the all-sector census, in sectors `4,5`, solve the Hermitian
definite pencil

```text
A=-W*[(V+V*)/2]W,
B=-W*[(M+M*)/2]W>0,
A e=lambda B e.
```

Require `15` negative and `10` positive generalized values with their frozen
sign margins.  Let `E_-` be the rank-`15` negative generalized space and
`P_L,P_-` the Euclidean projectors onto `L_H,E_-`.  For all `16` selected
cells record:

1. `||P_L-P_-||_2`;
2. `||L* A T_H||_2` and `||L* B T_H||_2`;
3. the inertia of `L* A L` and `T_H* A T_H`;
4. the leakage of the centered `Gamma,Omega` from `L_H`;
5. schedule/variant comparisons.

Equality is not inferred from dimension.  It requires projector equality,
zero cross blocks, negative definiteness on `L_H` and positive definiteness on
`T_H` under the frozen `10/100` bands.

## 6. Negative control

Construct a deterministic same-dimension control by replacing the first
Euclidean singular vector of `L_H` with the first vector of `T_H`, followed by
QR orthonormalization.  This control is fixed independently of the generalized
eigenvectors.  Its projector distance from `E_-` must be resolved nonzero; it
must not satisfy the complete longitudinal identity.

## 7. Outcomes and scope

Return

```text
NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_RESOLVED
```

only if all geometry, all-sector census, selected projector, cross-block,
sign, dynamic-invariance and negative-control gates pass.

Return

```text
NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_REFUTED
```

if any selected equality or invariance is resolved nonzero, and
`NEGATIVE_SHAPE_LONGITUDINAL_IDENTITY_OPEN` for unresolved margins.

A positive result is **DERIVED COMPUTATIONAL / STRUCTURAL**: the coarse
negative carrier is the action-weighted image of tangential vertex
displacements.  It is not an exact gauge theorem on the curved finite
background and does not authorize a quotient.  A negative result leaves its
scalar/vector/tensor identity **OPEN**.

Only the new targeted verifier and static guards may run.  The full suite is
excluded.
