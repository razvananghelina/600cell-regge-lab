# Prior-art gate: does 600-cell curvature holonomy kill the common flex seed?

Date: 2026-08-19

## Exact object and complete hypotheses

Use the regular piecewise-Euclidean boundary triangulation of the 600-cell:

```text
(f0,f1,f2,f3)=(120,720,1200,600),
five regular tetrahedra incident on every spatial edge.
```

At a homothetic Lorentzian product placement

```text
q_i=lambda p_i+tau n,   tau!=0,
```

each tetrahedral frustum has the accepted six-dimensional local
length-plus-strut flex kernel.  The accepted two-frustum theorem says that
shared-face compatibility contains no relative mode: after transporting
frames across the face, adjacent local flexes must be identical.

On the connected dual graph, a putative global flex is therefore determined
by one local seed.  Transport around a closed dual loop must return that seed
to itself.  This mission asks:

> What is the common fixed subspace of the affine Regge holonomies around
> two nonparallel edges of one base tetrahedron, restricted to the accepted
> six-dimensional local flex kernel?

This is an infinitesimal global-closure question at the regular homothetic
background.  It contains no Regge action, dust variation, Hessian spectrum
or continuum limit.

## Regge holonomy under test

A regular tetrahedron has interior dihedral angle

```text
theta=arccos(1/3).
```

Five tetrahedra meet at each 600-cell edge, so the spatial deficit is

```text
delta=2 pi-5 theta.
```

Parallel transport around that edge is the affine spatial rotation about
the edge line by `delta`, embedded trivially on the time direction.  The
proposed exact trigonometric controls are

```text
cos(delta)=241/243,
sin(delta)=22 sqrt(2)/243.
```

These values are disclosed predictions and must be rederived from the
regular tetrahedron Gram matrix and Chebyshev identities before any fixed
space is accepted.

For an affine isometry

```text
h(x)=L x+c,
```

the induced action on a Poincare Killing field `X(x)=A x+b` is

```text
A' = L A L^(-1),
b' = L b-A' c.
```

The flex seed closes around a hinge only if `(A',b')=(A,b)`.

## Disclosed fixed-space prediction

Choose two nonparallel edges sharing a vertex in one regular tetrahedron.
For one edge holonomy, the expected local fixed dimensions are

```text
static local kernel       2
    rotation about the edge + translation along it;

expanding local kernel    1
    the boost/translation combination along the edge.
```

For two nonparallel edge holonomies, the disclosed prediction is

```text
common fixed local kernel dimension 0
```

on both static and expanding strata.

As controls:

- zero deficit must leave all six local directions;
- the full ten-dimensional Poincare algebra should have a one-dimensional
  common fixed space for the two rotations, namely time translation;
- that time translation is absent from both accepted local kernels.

No holonomy fixed-space matrix has been evaluated while writing this gate.

## What the primary literature establishes

### KNOWN

- In Regge calculus, curvature is concentrated on codimension-two hinges and
  parallel transport around one hinge gives a simple rotation related to the
  deficit angle.  In three dimensions the intrinsic curvature/holonomy and
  deficit-angle relation is treated explicitly by Ariwahjoedi and Zen,
  [*(2+1) Regge Calculus: Discrete Curvatures, Bianchi Identity, and
  Gauss-Codazzi Equation*](https://arxiv.org/abs/1709.08373).
- Holonomies circling a single hinge are simple rotations in the orthogonal
  plane, and their angle relations encode discrete Bianchi identities:
  Ariwahjoedi and Zen,
  [*Contracted Bianchi Identity and Angle Relation on n-dimensional
  Simplicial Complex of Regge Calculus*](https://arxiv.org/abs/1807.11420).
- Coordinate/frame formulations of simplicial gravity define parallel
  transport and curvature from simplex metrics and face transitions:
  A. D'Adda,
  [*Simplicial Gravity with Coordinates*](https://arxiv.org/abs/2007.15361),
  and V. M. Khatsymovsky,
  [*Affine connection form of Regge calculus*](https://arxiv.org/abs/1509.04974).

These sources establish the holonomy framework.  They do not state the fixed
subspace of this repository's six frustum flexes.

### CONTROL

- exact 600-cell incidence: five tetrahedra per edge;
- exact regular-tetrahedron dihedral angle and nonzero deficit;
- the accepted one-frustum completeness theorem;
- the accepted two-frustum diagonal propagation theorem;
- origin covariance of the affine Poincare adjoint action.

### OPEN

- the exact fixed dimensions for one and two edge holonomies;
- whether two local hinge loops suffice to kill every global seed;
- whether the resulting zero seed is reproduced by a direct complete-carrier
  calculation rather than only by local holonomy generators;
- finite rather than infinitesimal uniqueness;
- differentiability and branch uniqueness of a globally reconstructed
  cellular action;
- physical Hessian, dynamics and continuum behavior.

No primary source found in the 2026-08-19 search performs this exact flex-
holonomy intersection.  External novelty is **OPEN**.

## Acceptance and kill boundary

- **Global infinitesimal rigidity candidate:** each one-edge control has its
  predicted nonzero fixed space, but two nonparallel edge holonomies have
  zero common fixed space on both strata.
- **Global underdetermination survives:** a positive common fixed space
  remains.
- **Structural failure:** the dimension depends on origin, edge orientation
  or an arbitrary development convention.
- **Control failure:** the deficit, affine adjoint, full-Poincare time-
  translation control or upstream theorems fail.

Even a zero common seed proves only that the previously identified local
flexes cannot extend to a global infinitesimal flex.  A mechanically
independent full-carrier audit is required before restoring a global
length-only anisotropic Hessian.
