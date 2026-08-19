# Protocol: exact affine-holonomy obstruction to a global frustum flex

Date: 2026-08-19

This protocol is committed before evaluating any holonomy fixed subspace.
It tests whether the one seed propagated by face gluing can close around two
nonparallel regular 600-cell edges.

## Frozen provenance

| input | SHA-256 |
|---|---|
| global-holonomy prior-art gate | `e5477823bc765d83cf812d393282ff8376c502d2967617226f42f1707474d056` |
| consolidated two-frustum result | `b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3` |
| primary two-frustum verifier | `52636ae59bd4e4568df175e32b7c3aeae4fbfbc3d475d255131b6db671c41ae7` |
| primary two-frustum artifact | `0e09c3f8f38c8158deff5b81bc6fe4d5d6dd685a24cce83e015fb95e3f26a70e` |
| adversarial two-frustum verifier | `b7a1f63e193aad50783929c8448ce99c18f1b50dc8e5ea27e3ed1102ec9dfa26` |
| adversarial two-frustum artifact | `0f8e70ef89b7fd5a8995349d40c77f6d3f637f2d9ce137ce2c9ff07b2fed2542` |
| complete 600-cell finite-Regge verifier | `4419c409c66488a246fbfdd2ff8ba265cf835884635738cfcf1bef8eab9ae5b2` |
| complete 600-cell finite-Regge artifact | `9f78212270f2dd2b3f73a2dd914f497a0914fc1f114c0483b1ea65b12036a025` |

The gluing artifacts must preserve their diagonal-only outcomes.  The
finite-Regge artifact must preserve `18/18` and `f=(120,720,1200,600)`; the
new verifier must independently check that five tetrahedra meet at each
edge rather than infer it from the f-vector alone.

## Exact base tetrahedron and hinge axes

Use the centered regular tetrahedron

```text
p0=( 1, 1, 1)   p1=( 1,-1,-1)
p2=(-1, 1,-1)   p3=(-1,-1, 1).
```

Choose the two nonparallel edges `(0,1)` and `(0,2)`, which share `p0`.
Their normalized directions may contain `sqrt(2)`, but the final rotation
matrices must simplify exactly over the rationals.

Derive the regular tetrahedron's interior dihedral cosine from its exact
face normals and require

```text
cos(theta)=1/3.
```

Using the exact fivefold incidence, derive

```text
delta=2 pi-5 theta,
cos(delta)=241/243,
sin(delta)=22 sqrt(2)/243 > 0.
```

The cosine must be obtained from the fifth Chebyshev polynomial and the sine
from its exact sign and `sin^2=1-cos^2`.  No floating angle or fitted deficit
is permitted.

For each edge, construct by Rodrigues' formula the spatial rotation `R_e`
through `delta` about its direction, then the affine rotation about the
actual edge line

```text
h_e(x)=R_e x+c_e,
c_e=p0-R_e p0.
```

Embed it in Poincare spacetime as `L_e=diag(R_e,1)` with zero time
translation.  Require exact orthogonality, determinant one, pointwise fixing
of both edge endpoints and noncommutation of the two linear rotations.

## Exact Poincare adjoint

In the disclosed basis of six Lorentz generators and four translations,
construct the ten-dimensional adjoint action

```text
A' = L_e A L_e^(-1),
b' = L_e b-A' c_e.
```

All coordinates and ranks must use exact SymPy algebra.  For each edge,
reversing loop orientation replaces the adjoint by its inverse and must leave
every fixed-space dimension unchanged.

Translate the entire base tetrahedron by `(2,-1,3)`, rebuild both affine
holonomies, and require the same fixed-space dimensions.  This is the
anti-development-convention control.

## Local flex kernels and frozen predictions

Use the analytic local parameter kernels at

```text
(lambda,tau)=(1,5),(2,5),(3,11).
```

For a zero-deficit identity holonomy, all six directions must remain fixed.

For either single nonzero edge holonomy require

```text
lambda=1:   fixed dimension 2,
lambda!=1:  fixed dimension 1.
```

For both nonparallel edge holonomies together require

```text
common fixed dimension 0
```

at all three representatives and under both loop orientations and the
origin shift.

On the full ten-dimensional Poincare algebra, the two-edge common fixed
space must instead have dimension one and equal exactly the time-translation
line.  This positive control shows that the holonomies do not simply have a
trivial common centralizer; the accepted local kernels specifically exclude
the surviving full-Poincare direction.

## Outcome hierarchy

1. `GLOBAL_FLEX_HOLONOMY_CONTROL_FAILED` if provenance, 600-cell incidence,
   dihedral/deficit derivation, affine-axis geometry, adjoint, orientation,
   origin or full-Poincare controls fail.
2. `GLOBAL_FLEX_SEED_KILLED_BY_HOLONOMY` if every one-edge prediction passes
   and the two-edge common local fixed space is zero on all strata.
3. `GLOBAL_FLEX_SEED_SURVIVES_HOLONOMY` if controls pass and any positive
   common local fixed space remains.
4. `GLOBAL_FLEX_HOLONOMY_OPEN` otherwise.

## Interpretation firewall

A zero seed establishes a **candidate exact infinitesimal global-rigidity
theorem** by combining local completeness, diagonal face propagation and
hinge closure.  It does not yet authorize a physical Hessian.  Acceptance
requires a mechanically independent audit on the complete 600-cell dual
complex, without replacing all loops by the two disclosed model
holonomies.

Even after that audit, finite branch uniqueness, a differentiable global
shape reconstruction, an action, dynamics and continuum convergence remain
separate questions.

Only the new verifier and static registry guards may be run.
