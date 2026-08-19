# Result: curvature holonomy closes the complete local frustum-flex family

Date: 2026-08-19

## Headline

The six infinitesimal length-preserving modes of one homothetic tetrahedral
frustum do **not** extend to a nonzero compatible mode on the complete regular
600-cell carrier.  Shared-face matching propagates one local seed, and the
actual Regge holonomies around the six edges of any base tetrahedron have zero
common fixed subspace on that seed.

The mechanically independent complete-dual-complex audit returned

```text
11/11 checks passed
ADVERSARIAL_GLOBAL_FLEX_SEED_KILLED
```

The exact status is:

```text
DERIVED EXACT, ADVERSARIALLY CORROBORATED:
the complete six-dimensional one-frustum flex family has no nonzero global
infinitesimal extension satisfying all fixed-frame face gluings on the regular
600-cell.
```

This is an infinitesimal global-closure result.  It is not yet finite branch
uniqueness, an action, an equation of motion, a graviton spectrum or a
continuum-limit claim.

## Complete theorem statement

Assume all of the following.

1. The spatial carrier is the regular Euclidean 600-cell boundary complex,
   with 120 vertices, 720 edges, 1200 triangular faces and 600 tetrahedra.
2. Each spacetime cell is a nondegenerate homothetic tetrahedral frustum with
   top vertices

   ```text
   q_i = lambda p_i + tau n,  tau != 0,
   ```

   and its infinitesimal variables are only the four top-vertex displacements.
3. In each cell all six top-edge squared lengths and all four strut squared
   lengths are preserved to first order.  By the previously derived relative
   Poincare theorem, this local kernel is exactly six-dimensional.
4. Adjacent cells identify their three shared top vertices in one fixed-frame
   geometry.  By the two-frustum theorem, two compatible local kernels agree
   diagonally; there is no independent relative face mode inside this family.
5. The propagated local Poincare seed is transported across faces by the
   unique vertex-matched Euclidean isometry, and it must close around every
   dual edge-star loop.

Then the only compatible global infinitesimal displacement in this local
kernel family is zero, modulo motions not contained in the accepted local
strut kernel.  Equivalently, the common fixed subspace of the six actual
base-edge holonomies, restricted to each accepted local kernel, has dimension
zero.

The qualification in the last sentence matters.  On the full ten-dimensional
Poincare algebra the same six holonomies retain exactly the one-dimensional
time-translation line.  The constraints of the homothetic frustum exclude
that positive-control direction.

## Provenance ledger

| stage | commit | result |
|---|---|---|
| prior-art gate | `18398e1` | discrete holonomy separated from the new closure claim |
| primary protocol | `19de8f0` | predictions frozen before model-loop ranks |
| primary registration | `1704b56` | exact Chebyshev/Rodrigues verifier frozen |
| primary artifact | `f1c28f2` | `12/12`, two nonparallel model loops kill the seed |
| adversarial protocol | `5ae2648` | complete actual dual-complex audit frozen |
| adversarial registration | `febe26c` | all 720 loops specified before evaluation |
| adversarial artifact | `5338374` | `11/11`, complete closure corroborated |

Frozen SHA-256 values:

```text
primary protocol
671cfcd02d902a8cc95969619c7ae9bdb3279efd4704ea210f00b0b337be66b1

primary verifier
9e4c13cf944283fbe473c318853ac951701abe6ac7147c78f525a1de071d7120

primary artifact
6852c4f0da3f747f178a697647bc0326a9668858ef414d0078668f2030875acf

adversarial protocol
9404da1f2a9ca5b1d7cf0038f81870aff8916a83f8014d7787ed14fa3915c325

adversarial verifier
54fa9775a2f14d708359167d3f8b81e03d985f24594b453f16028d9981d9be0d

adversarial artifact
f224fe123c882ccda97d4ca6ec67c9fd810d58ed8377c5afb457a1dec69f4b87
```

## Primary exact obstruction

The primary verifier derived the regular-tetrahedron dihedral data exactly,
without floating angles:

```text
cos(theta) = 1/3,
delta = 2 pi - 5 theta,
cos(delta) = 241/243,
sin(delta) = 22 sqrt(2)/243.
```

It constructed affine rotations about two disclosed nonparallel tetrahedron
edges and their exact adjoint action on the Poincare algebra.  For

```text
(lambda,tau) = (1,5), (2,5), (3,11),
```

one edge retained dimension two on the static stratum and dimension one on
each expanding stratum.  The two-edge common fixed dimension was zero in all
three cases.  Reversing the loops and shifting the affine origin changed no
decision.  The full-Poincare common fixed space was exactly time translation,
providing a nontrivial positive control.

## Mechanically independent complete audit

The adversarial verifier did not import the primary Rodrigues rotations.  It
independently rebuilt the exact golden-field carrier by clique incidence and
obtained

```text
f = (120,720,1200,600),
five tetrahedra per edge,
two tetrahedra per face,
connected 600-node dual graph.
```

It then constructed all 2400 directed face transitions by matching their
three shared vertices and reflecting the opposite apex.  Every transition
was an exact inverse isometry.

For every one of the 720 spatial edges it multiplied the five actual face
transitions around the edge star.  All 720 loops had

```text
det(R) = 1,
trace(R) = 725/243 = 1 + 2 cos(delta),
```

fixed both edge endpoints pointwise, and became their exact inverses when the
traversal was reversed.

For the six actual loops based at the lexicographically first tetrahedron the
fixed dimensions were

```text
(lambda,tau)=(1,5):  each single loop 2, all six 0
(lambda,tau)=(2,5):  each single loop 1, all six 0
(lambda,tau)=(3,11): each single loop 1, all six 0.
```

An odd relabelling of the base tetrahedron preserved the zero common fixed
space.  On the full Poincare algebra, all six loops again retained exactly
time translation.

## What changed conceptually

The earlier local no-go was too broad.  One isolated frustum really is
underdetermined by ten lengths, but those six local modes are not independent
variables on different cells.  Face matching propagates one seed, while
noncommuting curvature holonomies obstruct its global closure.

Thus both statements are simultaneously true:

```text
local cell:     six infinitesimal flexes;
complete 600-cell: no nonzero global flex in that complete local family.
```

This is the discrete analogue of the familiar fact that a parallel section
must be fixed by the holonomy group.  What is special here is the exact
600-cell calculation and its restriction to the analytically complete local
frustum kernel.

## Relation to prior art

Discrete face transports and curvature as products around hinges are standard
in connection formulations of Regge calculus; see, for example,
[Khatsymovsky](https://arxiv.org/abs/1509.04974).  Hessian zero modes and their
relation to discrete gauge symmetry are also established topics; see
[Dittrich--Freidel--Speziale](https://arxiv.org/abs/0707.4513).

The 600-cell has previously been used in Regge cosmology, including evolution
with more than one free variable by
[De Felice--Fabri](https://arxiv.org/abs/gr-qc/0106077).  The post-result search
found no primary source that states or computes the present all-720-loop
fixed-space theorem for the six local homothetic-frustum flexes.

That search absence is not proof of novelty.  External novelty remains
**OPEN** until a specialist literature review or peer review confirms it.

## Status ledger

| Claim | Status |
|---|---|
| One isolated homothetic tetrahedral frustum has six local length flexes | **DERIVED EXACT** |
| Those flexes are the restrictions of relative Poincare motions | **DERIVED EXACT** |
| Adjacent cells retain an independent relative flex inside that family | **REFUTED DERIVED EXACT** |
| The complete regular 600-cell has a nonzero global flex in that family | **REFUTED DERIVED EXACT, ADVERSARIALLY CORROBORATED** |
| All 720 actual hinge holonomies have the regular Regge deficit | **DERIVED EXACT** |
| The result is invariant under loop reversal and base relabelling | **DERIVED EXACT controls** |
| The complete nonlinear length reconstruction is locally unique | **OPEN** |
| A finite disconnected solution branch is absent | **OPEN** |
| A schedule-free anisotropic Regge action is now defined | **NOT YET** |
| The action Hessian has propagating physical modes | **NOT TESTED** |
| The continuum limit reproduces general relativity | **NOT TESTED** |
| The exact closure theorem is externally novel | **OPEN** |

## Next discriminating step

The holonomy result authorizes, but does not replace, an implicit global
reconstruction test.  Construct the complete constraint map from globally
identified vertex coordinates to the permitted length data and evaluate its
exact Jacobian modulo the expected global isometries.

- full transverse rank: an implicit local reconstruction map exists, so a
  schedule-free action and its anisotropic Hessian become legitimate;
- residual kernel: the surviving directions must be classified before any
  dynamics is claimed;
- rank only at the regular point but not nearby: the result may be a
  symmetry-stratum accident and does not authorize a generic evolution map.

Only after that regularity gate should the project return to the proposed
temporal-plus-spatial fluctuation equation and ask whether a physical wave
speed can be derived.
