# Adversarial protocol: complete dual-complex flex holonomy audit

Date: 2026-08-19

This protocol is committed after the primary two-model-holonomy result and
before reconstructing any face transition or dual loop.  The audit must not
import the primary Rodrigues rotations or replace actual dual loops by model
axes.

## Frozen inputs

| input | SHA-256 |
|---|---|
| primary holonomy protocol | `671cfcd02d902a8cc95969619c7ae9bdb3279efd4704ea210f00b0b337be66b1` |
| primary holonomy verifier | `9e4c13cf944283fbe473c318853ac951701abe6ac7147c78f525a1de071d7120` |
| primary holonomy artifact | `6852c4f0da3f747f178a697647bc0326a9668858ef414d0078668f2030875acf` |
| consolidated two-frustum result | `b5bb18c75ea1359d33b9985ad5816c21f437960c06f8c4eae793a3505509add3` |
| direct five-vertex gluing artifact | `0f8e70ef89b7fd5a8995349d40c77f6d3f637f2d9ce137ce2c9ff07b2fed2542` |

The primary artifact must retain `12/12` and
`GLOBAL_FLEX_SEED_KILLED_BY_HOLONOMY`.  The gluing inputs must retain their
diagonal-only conclusions.

## Independent complete carrier

Reconstruct the 120 exact golden-field vertices and derive all simplices by
clique incidence.  Require

```text
(f0,f1,f2,f3)=(120,720,1200,600),
two tetrahedra per face,
five tetrahedra per edge,
connected 600-node dual graph.
```

For every tetrahedron, assign canonical local coordinates by mapping its
sorted four global vertex labels to

```text
(1,1,1), (1,-1,-1), (-1,1,-1), (-1,-1,1).
```

This label-based choice is only a gauge.  All final loop invariants must be
independent of reversing the edge-star traversal.

## Face transitions constructed from vertex matching

For every directed pair of adjacent tetrahedra `T <- T'`:

1. keep the three shared global face vertices fixed in `T` coordinates;
2. reflect the apex of `T` across their affine face plane;
3. construct the unique affine map sending the canonical coordinates of
   `T'` to those four developed positions.

Require the linear part to be exactly orthogonal.  Do not insert a deficit
angle or axis rotation anywhere in this construction.

## All 720 actual hinge loops

For each spatial edge:

1. build the adjacency cycle of its five incident tetrahedra;
2. multiply the five directed face transitions around the cycle;
3. repeat with the reversed cycle.

Every resulting affine loop must:

- have orthogonal linear part and determinant one;
- fix both canonical edge endpoints pointwise;
- have exact linear trace

```text
1+2*(241/243)=725/243;
```

- give the same fixed spaces under reversed traversal.

The census must cover exactly all 720 edges.  This is the mechanically
independent derivation of the Regge deficit; the primary Chebyshev/Rodrigues
construction is used only as a frozen comparison afterward.

## Complete local-seed closure audit

Choose the lexicographically first base tetrahedron.  For each of its six
edges, use the actual loop based at that tetrahedron and construct its exact
ten-dimensional affine Poincare adjoint.

At `(lambda,tau)=(1,5),(2,5),(3,11)` require:

1. every single-edge fixed dimension agrees with the primary census
   (`2` static, `1` expanding);
2. the common fixed space of all six actual base-edge loops is zero;
3. the full ten-dimensional Poincare common fixed space is exactly the
   time-translation line;
4. reversing all six edge cycles changes no fixed dimension.

As a gauge attack, relabel the four vertices of the base tetrahedron by the
odd permutation `(0 1)` and rebuild its six loop matrices through the induced
canonical-coordinate change.  The common fixed dimensions must remain zero
on the local kernels and one on the full Poincare algebra.

## Outcome hierarchy

1. `ADVERSARIAL_GLOBAL_HOLONOMY_CONTROL_FAILED` if provenance, complete
   incidence, face-transition isometry, any of the 720 loop invariants,
   reversal or relabelling controls fail.
2. `ADVERSARIAL_GLOBAL_FLEX_SEED_KILLED` if the actual six base-edge loops
   have zero common fixed local seed on all three strata.
3. `ADVERSARIAL_GLOBAL_FLEX_SURVIVES` if controls pass but a positive common
   local fixed space remains.
4. `ADVERSARIAL_GLOBAL_HOLONOMY_OPEN` otherwise.

## Interpretation firewall

Corroboration proves an exact **infinitesimal global-rigidity theorem for the
specific local flex family**: no global flex can restrict to the complete
one-frustum kernels and satisfy all fixed-frame face gluings.

It still does not prove finite uniqueness or supply a practical global shape
reconstruction map.  A refined anisotropic action/Hessian is authorized only
after a separate implicit-reconstruction regularity test; it must not be
obtained by silently choosing one of the 24 staircase schedules.

Only this verifier and static registry guards may be run.
