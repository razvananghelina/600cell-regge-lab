# Preregistered protocol: direct-precision generalized-mode bundle

Date: 2026-08-18

Prior-art/framing commit: `6fdddc6`.

Status: **TARGET-DISCLOSED, PREREGISTERED BEFORE ANY DIRECT OLD/SHIFTED
GENERALIZED PROJECTOR, PROJECTOR DISTANCE, OR DIRECT `Gamma/Omega` LEAKAGE IS
COMPUTED.**

The preceding broad-error labels and approximate `10^-7` binary midpoint
diagnostics are known.  No characteristic root or spatial spectral target is
known or admissible here.

## 1. Frozen provenance

Require the following exact inputs and SHA-256 values:

| input | SHA-256 |
|---|---|
| first accepted tick JSON | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |
| second accepted tick JSON | `936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70` |
| third accepted tick JSON | `ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0` |
| geometry verifier | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| full-rank source | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| full tangent source | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| conformal source | `d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4` |
| `commons/cell600.py` | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |
| old centered JSON | `fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56` |
| old centered NPZ | `1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef` |
| shifted centered JSON | `265bd863de2365f19f7679373155fdaa23fb0bb3e75c221cfd9d9ec5b6ac2a47` |
| shifted centered NPZ | `c000f4fcae67e6c0648046878c2bd1ffd0616c38510ccf788c67cf99832397b8` |
| shifted direct-rank artifact | `86b53f228d6cfa7326a677d881463f1b849e76bc6c9ac2b0e8aa6fd427042944` |
| broad generalized-closure artifact | `53e046e2020a97fc992559546ce3d45479c0c0de7ce2e01322b09998ba85cf80` |
| shifted direct-rank source | `1b54cd25899037fc66c2b58e01ef3bac267c6ebf2c6917d2a05ac4ac0feed1c5` |
| broad generalized-closure source | `0a84c8ec4fab1c9626d5e4c711f89c6f9638cf37c15ef6b0050d6b66dfdde6c1` |

Require the upstream outcomes

```text
HOMOTHETIC_CANONICAL_LAPSE_SELECTED,
SECOND_HOMOTHETIC_TICK_ACCEPTED,
THIRD_HOMOTHETIC_TICK_ACCEPTED,
CENTERED_JACOBI_CERTIFIED,
SHIFTED_CENTERED_CERTIFIED,
SHIFTED_DIRECT_NEGATIVE_RANK_PERSISTS,
GENERALIZED_MODE_RECURRENCE_CLOSURE_CERTIFIED.
```

## 2. Complete direct reconstruction

Use `100` decimal digits for midpoint work and `80` decimal digits for Flint
balls.  Retain the frozen derivative steps

```text
operational_primary   1e-20,
operational_shadow    1e-15,
validation_primary    3e-20,
validation_shadow     3e-15.
```

For each schedule parity reconstruct directly from the action:

```text
slab 1: (0,a1,r1),
slab 2: (a1,a2,r2),
slab 3: (a2,a3,r3).
```

For sectors `4,5` and every derivative variant, project the three complete
Hessians, apply the already-authorized unique Hermitian projection
`(H+H*)/2`, build the tangent balls, and reconstruct the principal blocks
`S00,S01,S10,S11`.  Require every raw anti-Hermitian defect to lie inside the
complete four-variant Hessian family variation, every boundary twist to be
regular, and every principal-function recovery identity to contain zero
entrywise.

Form the two centered stencils without serialization:

```text
old:
  Kminus=S1,10,
  Kzero =S1,11+S2,00,
  Kplus =S2,01;

shifted:
  Kminus=S2,10,
  Kzero =S2,11+S3,00,
  Kplus =S3,01.
```

Then form

```text
M=(Kminus+Kplus)/2,
N=(Kplus-Kminus)/2,
V=Kminus+Kzero+Kplus,
Gamma=M^-1 N,
Omega=M^-1 V.
```

Every direct ball must overlap its same-cell broader serialized control.
Record raw interval-radius reduction factors for all five matrices.  Preserve
the upstream requirement that every direct `V` radius improves by more than
`100`; no new improvement threshold is imposed after seeing the other four.

## 3. Finite-difference family enclosure

A narrow fixed-step Flint radius controls arithmetic and source ball
propagation, not the `h -> 0` derivative limit.  For every centered time,
parity, sector, variant and matrix `X`, define

```text
epsilon_family(X_v)
  = max_w [ ||mid(X_v)-mid(X_w)||_2
            + radius_F(X_v)+radius_F(X_w) ]
    + arithmetic_floor.
```

Here `radius_F` is the Frobenius norm of the componentwise reenclosure radii.
This complete four-variant envelope is the admitted matrix error.  Do not use
the much smaller fixed-step interval radius by itself.

The statement remains conditional on this frozen finite-difference family;
it is not promoted to a formal interval theorem for the exact analytic
Hessian.

## 4. Shape carrier and generalized projector

Reconstruct the exact incidence carrier and deterministic sector bases.  Let
`U` span the rank-five conformal incidence image inside each 30-dimensional
minimal sector.  With the Hermitian midpoint `H_M`, obtain the
action-selected 25-dimensional shape carrier `W` as the nullspace of

```text
U* H_M.
```

Propagate incidence and row-nullspace errors exactly as in the broad
generalized-mode verifier, but replace every serialized source error by the
family envelope above.

On the shape carrier form

```text
A=-W* H_V W,
B=-W* H_M W.
```

Require `B` positive-definite-resolved and the generalized negative/positive
cluster gap separated after its complete error.  Solve only through the
standard Hermitian-definite driver and lift the first 15 generalized
eigenvectors.  QR is allowed only to obtain an orthonormal basis for their
unchanged span.  Record the Hermitian projector `P` and propagate the same
conservative pencil/projector bound used upstream.

Required census:

```text
32 projectors = 2 times * 2 parities * 2 sectors * 4 variants.
```

Each must retain exactly `15 negative + 10 positive` generalized eigenvalues.

## 5. Identity and closure tests

### Projector identity

For each of the 16 parity/sector/variant pairs, compare

```text
d=||P_shifted-P_old||_2,
epsilon_d=eta_old+eta_shifted+arithmetic_floor.
```

Classify with the frozen bands:

```text
d <= 10 epsilon_d     GENERALIZED_COMMON_FIBER_RESOLVED,
d > 100 epsilon_d     GENERALIZED_ROTATED_FIBER_RESOLVED,
otherwise             GENERALIZED_FIBER_IDENTITY_OPEN.
```

No eigenvector matching or fitted alignment is permitted.

### Local and cross-slice closure

For every direct projector `P_t`, test both direct operators at both times:

```text
R(P_t,X_s)=(I-P_t) X_s P_t,
X in {Gamma,Omega},
t,s in {old,shifted}.
```

This gives

```text
128 leakages,
64 local  (t=s),
64 cross  (t!=s).
```

Use

```text
epsilon_R = epsilon_X
          + (2 eta_P+eta_P^2)(||X||+epsilon_X)
          + arithmetic_floor
```

and classify

```text
||R|| <= 10 epsilon_R     LEAKAGE_ZERO_CONSISTENT,
||R|| > 100 epsilon_R     LEAKAGE_NONZERO_RESOLVED,
otherwise                 LEAKAGE_OPEN.
```

Local closure asks whether the action-selected fiber is a recurrence fiber at
its own slice.  Cross closure asks the stronger fixed-carrier question.  Keep
the two ledgers separate.

## 6. Frozen outcome hierarchy

Apply the first matching outcome:

1. `DIRECT_GENERALIZED_BUNDLE_CONTROL_FAILED` on any provenance, branch,
   Hermitian-family, twist, principal-identity, overlap, precision, carrier,
   positivity, gap, finiteness or census failure;
2. `DIRECT_GENERALIZED_LOCAL_CLOSURE_REFUTED` if any local leakage is
   nonzero-resolved;
3. `DIRECT_GENERALIZED_LOCAL_CLOSURE_OPEN` if none is nonzero-resolved but any
   local leakage is open;
4. `DIRECT_GENERALIZED_BUNDLE_NONUNIFORM` if the 16 projector-identity labels
   are not all identical;
5. `DIRECT_GENERALIZED_BUNDLE_ROTATION_RESOLVED` if all 16 identity labels are
   `GENERALIZED_ROTATED_FIBER_RESOLVED`;
6. `DIRECT_GENERALIZED_BUNDLE_IDENTITY_OPEN` if all 16 identity labels are
   `GENERALIZED_FIBER_IDENTITY_OPEN`;
7. `DIRECT_GENERALIZED_COMMON_CROSS_CLOSURE_REFUTED` if all projectors are
   common but any cross leakage is nonzero-resolved;
8. `DIRECT_GENERALIZED_COMMON_CROSS_CLOSURE_OPEN` if all projectors are common,
   no cross leakage is nonzero-resolved, but at least one is open;
9. `DIRECT_GENERALIZED_COMMON_BUNDLE_RESOLVED` only if all 16 projectors are
   common and all 128 local/cross leakages are zero-consistent.

Outcome 5 is a clean result, not a failure: it says that the dynamical mode
space is a genuine bundle rather than one fixed subspace.  It does **not**
authorize the polar/direct rotation as physical transport.  A subsequent
protocol would have to compare any mathematical transport with the
action-generated tangent before using it.

No reduced companion, root count, Lyapunov exponent, dispersion, inertia,
mass or limiting-speed claim is allowed in this verifier.

## 7. Transparent post-first-result diagnostic amendment

The first two executions were byte-identical, returned `13/13`, and produced
the artifact preserved in commit `6901ebd` with SHA-256
`a7f6f915b9284905ad1931131edaa5cd2402dd3b13d1161be12e4201252641a7`.
The mechanical outcome was
`DIRECT_GENERALIZED_COMMON_BUNDLE_RESOLVED`.

Before choosing the next mission, add only the following diagnostic fields to
the artifact:

- the complete family error for each direct matrix;
- `eta_k`, the shape-row error and `eta_s`;
- the restricted `M,V` errors and the positive kinetic lower bound;
- the three additive contributions to `eta_P`: shape, generalized-eigenspace,
  and kinetic-metric terms.

This amendment changes no matrix, projector, residual, threshold, label,
outcome branch or check.  Its sole purpose is to identify which preregistered
error term limits the already disclosed common-bundle classification.  The
old artifact remains the immutable first-result record.
