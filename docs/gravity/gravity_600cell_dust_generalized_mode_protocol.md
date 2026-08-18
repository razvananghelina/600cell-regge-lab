# Preregistered protocol: generalized negative-mode recurrence closure

Date: 2026-08-18

Prior-art/framing commit: `52ec90a`.

Status: **TARGET-DISCLOSED, PREREGISTERED BEFORE ANY GENERALIZED PROJECTOR OR
`Gamma/Omega` LEAKAGE IS COMPUTED.**

## Frozen inputs

Use the same exact-provenance old/shifted centered archives, direct shifted
rank artifact, geometry, conformal helper, sector basis and commons source
frozen in commit `07de221`.  Their hashes remain those listed there.  Add the
negative-fiber transport artifact

```text
gravity_600cell_dust_negative_fiber_transport.json
d630bf07066f88c35eee5a62a80ec1f43399a95ea882a43528289220c67f4599.
```

Require its outcome `NEGATIVE_FIBER_TANGENT_CLOSURE_REFUTED` and the exact
`A,C` zero / `B,D` nonzero pattern.  No numerical projector from that artifact
is reused; all generalized fibers are rebuilt.

## Complete census

Use old and shifted centered slices, both parities, sectors `4,5` and all four
derivative variants:

```text
32 generalized projectors,
16 old/shifted identity comparisons,
64 recurrence leakages = 32 projectors * {Gamma,Omega}.
```

## Pencil and projector construction

Rebuild the rank-`25` action-selected shape carrier.  On it form the Hermitian
definite pencil

```text
A=-V_S,   B=-M_S.
```

Require `B` positive-definite-resolved.  Solve `A v=lambda B v` with the
standard Cholesky/Hermitian definite driver.  Require exactly `15` negative
and `10` positive midpoint generalized eigenvalues and a separated gap

```text
gap=lambda_15-lambda_14.
```

With inherited restricted errors `epsilon_A,epsilon_B`, define

```text
b_lower = lambda_min(B)-epsilon_B,
epsilon_pencil = epsilon_A/b_lower
                + ||A|| epsilon_B/(lambda_min(B) b_lower)
                + arithmetic_floor,

eta_eig = 2 epsilon_pencil/(gap-2 epsilon_pencil)
          + arithmetic_floor,

eta_P = 2 eta_shape
        + sqrt(lambda_max(B)/b_lower) eta_eig
        + epsilon_B/b_lower
        + arithmetic_floor.
```

Require all denominators positive.  Orthonormalize the lifted negative
generalized eigenvectors in the ambient Euclidean coordinate metric only to
form their unique Hermitian span projector.  This QR changes no subspace.

The formula is a deliberately conservative propagated bound conditional on
the frozen derivative family; it is not advertised as an optimal generalized
Davis--Kahan theorem.

## Identity and recurrence tests

Compare old/shifted generalized projectors under the identity seam map with
the existing `10/100` common/rotated/open classification.

For each old and shifted cell, load the corresponding normalized centered
operators and compute

```text
R_Gamma=(I-P) Gamma P,
R_Omega=(I-P) Omega P.
```

For `X` equal to `Gamma` or `Omega`, use

```text
epsilon_R = epsilon_X
          + (2 eta_P + eta_P^2)(||X||+epsilon_X)
          + arithmetic_floor.
```

Classify leakage with the same frozen thresholds:

```text
ZERO_CONSISTENT  ||R|| <= 10 epsilon_R,
NONZERO_RESOLVED ||R|| > 100 epsilon_R,
OPEN             otherwise.
```

## Frozen outcome hierarchy

1. `GENERALIZED_MODE_CONTROL_FAILED` on any provenance, carrier, kinetic,
   gap, error, finiteness or census failure;
2. `GENERALIZED_MODE_RECURRENCE_CLOSURE_REFUTED` if any of the `64` leakages
   is nonzero-resolved;
3. `GENERALIZED_MODE_RECURRENCE_CLOSURE_OPEN` if none is nonzero-resolved but
   any is open;
4. `GENERALIZED_MODE_RECURRENCE_CLOSURE_CERTIFIED` only if all `64` leakages
   are zero-consistent.

Report identity-fiber labels and `Gamma/Omega` labels separately.  No fitted
alignment, momentum graph, omitted operator, physical mass or full-suite run.
