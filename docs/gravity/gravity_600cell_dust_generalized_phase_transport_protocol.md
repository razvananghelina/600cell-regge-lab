# Protocol: canonical phase transport of the generalized fibers

Date: 2026-08-18

Prior-art/framing gate commit: `60aafe1`.

## 1. Frozen question

For each of the two disclosed symmetry sectors, both schedule parities and all
four derivative schedules, test whether the exact second-slab canonical
tangent transports the full cotangent lift of the old rank-15 generalized
fiber into the corresponding shifted lift.

For projectors `P0,P1` and tangent `T_2=[A B;C D]`, the fixed phase projectors
are

```text
Q0 = diag(P0, conjugate(P0)),
Q1 = diag(P1, conjugate(P1)).
```

No alternative graph, basis alignment or optimization may be introduced.

## 2. Frozen inputs and reconstruction

Require the residual-certified verifier and artifact with hashes

```text
verifier  ccf2ebe03c6e39c3d6e6b538d1c02d278804553987d65db0eeb67fce7936ca5a
artifact  3244185127aecf7c9a44261cced0be521c9dc42bf8e44f909d8a0ce10a96eadf
```

and outcome `RESIDUAL_FINITE_FAMILY_ROTATION_RESOLVED`.  Replay it completely
to recover its high-precision projectors and require `10/10` with the same
artifact hash.

Independently reconstruct only the physical second slab `(a1,a2,r2)` for both
parities at `100` decimal digits with Flint ball arithmetic at `80` decimal
digits.  For sectors `4,5` and all four schedules, rebuild the Hermitian action
Hessian, eliminate the 35 internal orbit variables and construct the exact
canonical tangent ball with the unique boundary map.  Require:

- all branch and principal-function identities;
- nonzero internal determinant;
- the identity boundary map;
- all entries of `T_2* Omega T_2-Omega` to contain zero;
- exactly 16 finite `60 x 60` tangent balls.

Do not load the old binary tangent archive as the scientific matrix.

## 3. Complete leakage census

For every matched `(parity,sector,schedule)` cell compute

```text
R_A = (I-P1)             A P0,
R_B = (I-P1)             B conjugate(P0),
R_C = (I-conjugate(P1))  C P0,
R_D = (I-conjugate(P1))  D conjugate(P0),

R_full = (I-Q1) T_2 Q0.
```

This gives exactly `64 = 16*4` block residuals and `16` full phase residuals.

For every tangent block `X`, let `epsilon_X` be twice the Frobenius norm of
the Flint component radii plus

```text
1e-70 * max(1,||X||).
```

With the already certified projector errors `eta_0,eta_1`, use

```text
epsilon_R = epsilon_X
          + (eta_0 + eta_1 + eta_0 eta_1)
            (||X|| + epsilon_X)
          + floor.
```

Use the same formula on the full tangent and phase projectors.  Record
midpoint norm, error, error units and label for every residual.

As a construction control, require the direct full residual matrix to equal
the `2 x 2` matrix assembled from `R_A,R_B,R_C,R_D` to relative `1e-60` at the
high-precision midpoint.  Require the full residual label to be compatible
with its four block labels: a nonzero-resolved block forces a
nonzero-resolved full residual, while four zero-consistent blocks force a
zero-consistent full residual; otherwise the cell is open.

## 4. Frozen labels

For residual norm `r` and complete error `e`, assign only

```text
r <= 10 e       LEAKAGE_ZERO_CONSISTENT,
r >  100 e      LEAKAGE_NONZERO_RESOLVED,
otherwise       LEAKAGE_OPEN.
```

No block or schedule may be dropped.

## 5. Frozen outcome hierarchy

Use the first applicable branch:

1. any provenance, reconstruction, symplecticity, rank, finiteness, assembly,
   label-compatibility or census control fails:
   `GENERALIZED_PHASE_TRANSPORT_CONTROL_FAILED`;
2. any of the 64 block or 16 full residuals is nonzero-resolved:
   `GENERALIZED_PHASE_TRANSPORT_REFUTED`;
3. none is nonzero-resolved but any is open:
   `GENERALIZED_PHASE_TRANSPORT_OPEN`;
4. all 80 are zero-consistent:
   `GENERALIZED_PHASE_TRANSPORT_CERTIFIED`.

A refutation closes only the unrestricted cotangent lift.  It does not close
the separately preregistered possibility of a canonical transported
intersection or Lagrangian graph.

## 6. Deliverable and exclusions

Write one deterministic JSON artifact with all hashes, reconstruction
controls, 64 block records, 16 full records, counts, outcome and status
ledger.  Register the verifier before its first scientific execution and run
it twice with a byte-identical artifact.

Run no unrelated verifier and no full suite.  Compute no intersection, graph,
Riccati solution, reduced root, dispersion, mass, limiting speed, graviton or
particle-inertia quantity in this mission.
