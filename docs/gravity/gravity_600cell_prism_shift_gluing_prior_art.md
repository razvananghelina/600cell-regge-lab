# Prior-art gate: global gluing of the static prism shift modes

Date: 2026-08-19

Status: **completed before protocol and before the carrier census**.

## 1. Exact object and complete hypotheses

Let `K` be either the boundary complex of the regular 600-cell or the already
certified projected rank-edgewise refinement of that complex.  For every
spatial tetrahedron `sigma` form one equal-scale Lorentzian polytopal cell

```text
sigma x I = Delta_3 x I.
```

The bottom and top copies of `K` have identical intrinsic spatial edge
lengths.  Every vertical edge has the same fixed timelike squared length
`-rho`, with `rho>0`.  No staircase diagonal, vertex order, embedding
coordinate, normal, shift or cross-edge length is supplied as input.

The local rigidity mission established that, after the bottom tetrahedron is
pinned, one such equal-scale cell has exactly three length-preserving modes:
the top tetrahedron is translated by a spatial tangential vector while the
timelike norm of the translation is fixed.  The present question is whether
those local modes survive the following **full three-face shape-matching
hypothesis**:

> If two spatial tetrahedra share a triangle, the complete induced metric of
> their common lateral triangular prism must agree under the labelled face
> identification.

The calculation is infinitesimal in the spatial shift around a nonzero
timelike translation, but the anticipated construction also has a finite
version.  It does not impose the Regge equations of motion, classify a mode
as gauge, or assume a continuum ADM interpretation.

## 2. The framing trap

A first tempting argument was:

```text
all four top vertices of one cell move together;
adjacent cells share three top vertices;
therefore adjacent translations are equal;
the dual graph is connected;
therefore only one global translation survives.
```

This is wrong for an intrinsic polytopal complex.  Adjacent four-cells do not
come with one common ambient affine frame.  Their shared triangular prism
sees only the restriction of the shift covector to the two-dimensional
triangle.  The components normal to that triangle inside the two spatial
tetrahedra need not be equal.

This correction is made before the protocol.  The proposed calculation must
therefore impose equality of induced face Gram matrices, not equality of
ambient vectors.

## 3. Local face metric

Choose two independent edge vectors `e1,e2` of a shared bottom triangle and
write the local cell translation as

```text
t = N n + s,
```

where `n` is timelike normal to the spatial tetrahedron and `s` is its
spatial tangential shift.  Relative to one bottom vertex, the induced Gram
matrix of the six-vertex lateral triangular prism is determined by

```text
e_i . e_j,   e_i . s,   t . t = -rho.
```

The first entries are fixed by the spatial slice and the last is fixed by the
vertical edge length.  Hence labelled shape matching across the face is
equivalent, at first order, to equality of the two numbers

```text
e1 . s,  e2 . s.
```

After metric duality, the local shifts are therefore constant one-forms on
tetrahedra whose tangential traces agree across triangular faces.

## 4. Mathematical prior art

The relevant algebra is not new.  Piecewise polynomial differential forms
with trace continuity form discrete de Rham complexes.  In particular, the
gradient of a continuous piecewise-affine scalar is a constant one-form on
each tetrahedron, and its tangential trace agrees across every face.  Closed
discrete one-forms modulo gradients represent first cohomology.

The standard finite-element exterior-calculus references are:

- D. N. Arnold, R. S. Falk and R. Winther, *Finite element exterior
  calculus, homological techniques, and applications*, Acta Numerica 15
  (2006), 1--155, DOI `10.1017/S0962492906210018`;
- D. N. Arnold, R. S. Falk and R. Winther, *Finite element exterior
  calculus: from Hodge theory to numerical stability*, Bull. Amer. Math.
  Soc. 47 (2010), 281--354, arXiv:`0906.4325`, DOI
  `10.1090/S0273-0979-10-01278-4`, especially the discrete de Rham and
  Whitney-complex discussion in Sections 5.1--5.6.

The present lowest-degree subspace can also be described without finite
element terminology.  Integrating a constant local one-form along each edge
gives a global edge cochain.  Face trace matching makes the edge value
independent of the incident tetrahedron, and local exactness makes its sum
around every triangle zero.  Thus it lies in

```text
ker(d1 : C^1(K;R) -> C^2(K;R)).
```

Conversely every such cocycle reconstructs one constant one-form per
nondegenerate tetrahedron.  Therefore the expected dimension is

```text
rank(d0) + b1 = V - b0 + b1.
```

For connected triangulated `S3`, `b0=1,b1=0`, giving `V-1`.  This is an
analytic prediction to be falsified by the registered matrix calculation,
not a number inferred from its output.

## 5. Regge-calculus boundary

Exact lapse/shift gauge generators are known on flat linearized Regge
backgrounds.  Hoehn shows that a `1-4` move creates four lapse/shift
variables and their vertex-displacement generators:

- P. A. Hoehn, *Canonical linearized Regge Calculus: counting lattice
  gravitons with Pachner moves*, Phys. Rev. D 91, 124034 (2015),
  arXiv:`1411.5672`, DOI `10.1103/PhysRevD.91.124034`.

That interpretation cannot be imported here.  Curvature generically breaks
exact vertex-displacement symmetry and turns constraints into
pseudo-constraints:

- B. Bahr and B. Dittrich, *(Broken) Gauge Symmetries and Constraints in
  Regge Calculus*, Class. Quantum Grav. 26, 225011 (2009),
  arXiv:`0905.1670`, DOI `10.1088/0264-9381/26/22/225011`.

Likewise, matching only areas or other incomplete face data can leave
non-metric shape mismatch.  Full shape matching is load-bearing here:

- S. K. Asante, B. Dittrich and H. M. Haggard, *The Degrees of Freedom of
  Area Regge Calculus: Dynamics, Non-metricity, and Broken
  Diffeomorphisms*, Class. Quantum Grav. 35, 135009 (2018),
  arXiv:`1802.09551`, DOI `10.1088/1361-6382/aac58d`.

## 6. KNOWN / CONTROL / OPEN

### KNOWN

- gradients of continuous `P1` potentials are piecewise-constant closed
  one-forms with matching tangential traces;
- discrete first cohomology measures closed forms not obtained this way;
- exact Regge lapse/shift gauge symmetry is background-dependent;
- incomplete face data can fail shape matching.

### CONTROL

- one equal-scale tetrahedral prism has exactly three local tangential modes;
- the 600-cell has `f=(120,720,1200,600)` and Betti numbers `(1,0,0,1)`;
- the projected rank-edgewise carrier is a subdivision of the same `S3` and
  has `f=(19680,134880,230400,115200)`;
- a single tetrahedron must give a three-dimensional local space;
- a complex with nonzero `b1` must supply a negative control in which the
  `V-1` formula fails by exactly `b1`.

### OPEN

- whether the literal face-Gram matching matrix on the 600-cell has dimension
  `119` and coincides exactly with the image of the vertex coboundary;
- whether the finite fixed-strut construction integrates every sufficiently
  small vertex potential;
- the number of invariant modes under the declared spatial symmetry action;
- whether the surviving modes are gauge, constrained or physical after a
  Regge action and its lapse/shift equations are included.

## 7. Decision boundary

If the full face-matching kernel is larger or smaller than the cocycle space,
the proposed reduction is wrong and the discrepancy is the result.

If the kernel equals `im(d0)` on the 600-cell, then the local three modes do
survive globally, but only as a longitudinal scalar-potential sector of
dimension `119`.  This would be a **DERIVED KINEMATIC** statement, not a
proof of ADM shift or gauge symmetry.

If finite cells reconstructed from two different nonconstant potentials have
the same declared bottom, top and vertical natural edge lengths while their
four-volumes or lateral face Gram data differ, then the static schedule-free
length-only cellular geometry is globally underdetermined.  That is a clean
negative for a canonical static Hessian from those lengths alone.

No action spectrum, desired dimension or physical target may be examined in
this mission.

