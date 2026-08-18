# Protocol: exact rank of the transported negative-stiffness intersection

Date: 2026-08-18

Prior-art/framing gate commit: `65419f6`.

Status: **TARGET-DISCLOSED, PREREGISTERED BEFORE ANY NEW INTERSECTION
SPECTRUM IS COMPUTED.**

## 1. Frozen inputs and replay

Require these SHA-256 hashes:

```text
high-precision phase source
  9c4c36b463a8faaa8d40b7db1b6b1852e3c04155c1b6ada4d02fbda747f6fcf3
high-precision phase artifact
  45eb9a3e80ead758d9b3c2f8e1eccff44b06e2759251ab00c447aa53e6705743
adversarial phase artifact
  c33615ac6d0f3133e53077f46c5ee766b9c633d4d64c32124c24839c9c84c880
binary negative-fiber source
  f462e507500d7f02ecf799f0d4b320e05795216a36a0d10eb908d6dc67b48181
binary negative-fiber artifact
  d630bf07066f88c35eee5a62a80ec1f43399a95ea882a43528289220c67f4599
```

Replay the accepted phase verifier in-process and require `7/7`, outcome
`GENERALIZED_PHASE_TRANSPORT_REFUTED`, the unchanged artifact hash, `32`
exact source-ball cells and `16` exact tangent balls.  Generalized projectors
and generalized transported-intersection singular values are controls only
and are not used to construct the negative projectors.

Replay the accepted binary negative-fiber verifier and require `8/8`, outcome
`NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED`, and `32` committed binary negative
projectors.  They are overlap controls only, never decisive inputs.

Use `mpmath` at `100` decimal digits, the already reconstructed Flint source
and tangent balls at `80` decimal digits, and a frozen arithmetic floor
`1e-70`.  No binary tangent archive may provide the decisive tangent.

## 2. High-precision negative projectors

For every time (`old,shifted`), parity, sector (`4,5`) and all four derivative
schedules, reconstruct the rank-`5` conformal carrier `U`, the exact midpoint
source matrices `M,V` and their complete ball-radius Frobenius bounds.

Rebuild the `25`-dimensional shape carrier as

```text
S = kernel(U* M)
```

with the frozen lexicographic pivot convention and Euclidean orthonormal
basis `W`.  If `g_S` is the smallest singular value of `U* M`, define the
shape-projector error

```text
eta_S = 2 epsilon_M/(g_S-2 epsilon_M) + floor
```

only when `g_S>2 epsilon_M`.

On this carrier form the ordinary, not generalized, Hermitian stiffness

```text
A = -W* V W
```

and propagate

```text
epsilon_A = epsilon_V
          + 2 eta_S (||V|| + epsilon_V)
          + floor max(1,||V||).
```

Diagonalize `A` at high precision.  Require exactly `15` midpoint-negative and
`10` midpoint-positive eigenvalues, both middle signs resolved against
`100 epsilon_A`, and the cluster gap

```text
g_A = lambda_15-lambda_14 > 2 epsilon_A.
```

For midpoint negative eigenbasis `Z_-` and complement `Z_+`, include both the
numerical cross-residual and source uncertainty:

```text
r_eig = ||Z_+* A Z_-||,
eta_res = sin(0.5 atan(2 r_eig/g_A)),
eta_src = 2 epsilon_A/(g_A-2 epsilon_A) + floor,
P^- = W Z_- Z_-* W*,
eta_P = 2 eta_S + eta_res + eta_src + floor.
```

Require Hermiticity, idempotence, shape-null and orthogonality controls.  The
high-precision `P^-` must overlap the matching committed binary projector
within that binary projector's conservative `eta`.  This overlap does not
replace the new error bound.

The carrier Hilbert metric selecting this spectral projector is frozen.  No
alternative metric or basis optimization is searched.

## 3. Complete blind intersection census

For every one of the `16` matched cells set

```text
Q_0 = diag(P^-_old, conjugate(P^-_old)),
Q_1 = diag(P^-_shifted, conjugate(P^-_shifted)),
R^- = (I-Q_1) T_2 Q_0.
```

For tangent midpoint norm `||T_2||`, full tangent ball error `epsilon_T` and
projector errors `eta_0,eta_1`, use the complete residual bound

```text
epsilon_R = epsilon_T
          + (eta_0+eta_1+eta_0 eta_1)(||T_2||+epsilon_T)
          + floor max(1,||T_2||).
```

Compute all `60` midpoint singular values at `100` digits and assign only

```text
s <= 10 epsilon_R   SINGULAR_ZERO_CONSISTENT
s > 100 epsilon_R   SINGULAR_NONZERO_RESOLVED
otherwise           SINGULAR_OPEN.
```

Record all values, bounds, error units and labels.  No cell, schedule or
singular value may be dropped.

## 4. Structural rank theorem and controls

Since the exact operator has right factor `Q_0` of exact rank `30`,

```text
rank(R^-) <= 30.
```

Require the lower `30` midpoint singular values to be zero-consistent and no
cell to contain more than `30` nonzero-resolved singular values.  If `30`
values are nonzero-resolved, Weyl's inequality plus the structural upper
bound certifies

```text
rank(R^-)=30,
dim K^-=30-rank(R^-)=0.
```

Fewer than `30` resolved singular values do **not** certify a positive
intersection; that cell remains open absent a separately preregistered exact
upper-rank proof.

Before the scientific census, require two exact synthetic controls in the
same `60`-dimensional layout:

1. `Q_0=Q_1=diag(I_30,0)`, `T=I`: leakage rank `0`, intersection dimension
   `30`;
2. `Q_0=diag(I_30,0)`, `Q_1=diag(0,I_30)`, `T=I`: leakage rank `30`,
   intersection dimension `0`.

These controls test the rank/nullity interpretation, not the physical data.

## 5. Frozen outcome hierarchy

Use the first applicable branch:

1. provenance, replay, reconstruction, sign/gap resolution, projector
   overlap, finiteness, ordering, controls, census or structural bound fails:
   `NEGATIVE_TRANSPORTED_INTERSECTION_CONTROL_FAILED`;
2. all `16` cells have exactly `30` nonzero-resolved singular values:
   `NEGATIVE_TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL`;
3. otherwise:
   `NEGATIVE_TRANSPORTED_INTERSECTION_DIMENSION_OPEN`.

Branch 2 closes the Hilbert-metric negative-fiber phase route.  Branch 3
authorizes no graph, Lagrangian or physical-mode claim.

## 6. Adversarial acceptance gate

This primary verifier cannot by itself consolidate a result.  After its
accepted execution, preregister and run a mechanically independent
replication which does not reuse the decisive high-precision projector/SVD
implementation.  It must include known full- and zero-intersection controls,
basis-gauge stress and a legitimate convention/precision stress.  Any
disagreement leaves the result **OPEN**.

## 7. Deliverable and exclusions

Write a deterministic JSON artifact containing hashes, all `32` projector
records, `16` cell summaries, all `960` singular-value records, controls,
counts, outcome and a status ledger.  Register the verifier before its first
scientific execution and run it twice with a byte-identical artifact.

Run no full suite and no unrelated verifier.  Compute no optimized alignment,
graph, symplectic restriction, propagator, root, dispersion, graviton, mass,
inertia or limiting speed in this mission.
