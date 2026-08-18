# Prior-art gate: canonical conformal subspace versus the Regge kinetic form

Date: 2026-08-18

## Exact object, carrier and hypotheses

Let `K` be the literal spatial boundary triangulation of the regular
600-cell used by the committed dust slabs.  It has

```text
|V(K)| = 120,   |E(K)| = 720,
```

and contains triangular faces.  The position carrier is the real vector
space of logarithmic squared spatial edge magnitudes

```text
Q = R^E,   q_uv = delta log |ell_uv^2|.
```

Define, without coefficients to choose, the unsigned vertex--edge incidence
map

```text
C : R^V -> Q,
(C sigma)_uv = sigma_u + sigma_v.
```

The irrelevant overall factor depends on whether the vertex variable scales
length, squared length or the metric.  It does not change `im C`, its rank,
or any restricted-sign conclusion.

Let

```text
H_p = (M_p + M_p*)/2
```

be the operational Hermitian centered coefficient already committed for
each staircase schedule `p=even,odd`.  In the fixed logarithmic edge
coordinates, the blind census found full-carrier inertia

```text
(n_+(H_p), n_-(H_p), n_0(H_p)) = (120, 600, 0).
```

The new primary question is not whether two dimensions happen to agree.  It
is whether the restriction

```text
G_p = C^T H_p C
```

is positive definite.  If `rank C=120`, positive definiteness makes
`im C` a maximal positive subspace because the complete positive index is
exactly `120`.  Reversing the overall action sign reverses both words
"positive" and "negative" but leaves the statement "the conformal image
carries the unique minority inertia count" unchanged.

All statements are on the fixed carrier, fixed dust background, fixed
literal identification of time fibres, two schedules, seven minimal binary-
tetrahedral sectors and four derivative variants.  They are not refinement
or continuum statements.

## Primary prior art

- [Glickenstein, *Discrete conformal variations and scalar curvature on
  piecewise flat two and three dimensional manifolds*](https://arxiv.org/abs/0906.1560)
  defines a perpendicular-bisector conformal structure by
  `ell_ij=exp(u_i+u_j)L_ij`.  Its tangent is
  `delta log ell_ij = delta u_i + delta u_j`; hence in logarithmic squared
  lengths it gives exactly the image of `C`, up to an overall factor.
  The paper treats piecewise-flat three-manifolds and variations of the
  Einstein--Hilbert--Regge functional.  Thus the vertex-scaling map itself is
  **KNOWN**, not a discovery of this project.
- [Hartle--Miller--Williams, *Signature of the Simplicial
  Supermetric*](https://arxiv.org/abs/gr-qc/9609028) derive and study the
  Lund--Regge metric on squared-edge configuration space.  They prove at
  least one physical timelike direction but also exhibit degeneracy,
  signature change and additional physical negative directions on other
  triangulations.  Generic Regge calculus therefore does not force the
  present `120:600` result.
- [Williams, *Recent Progress in Regge
  Calculus*](https://arxiv.org/abs/gr-qc/9702006) reviews the same warning:
  simplicial supermetric signature and gauge structure depend on the
  triangulation and background.

No located primary source computes the centered two-slab dust coefficient
on the fixed 600-cell, restricts it to this `120`-dimensional vertex-conformal
image, or establishes the result of the test below.  External novelty is
**OPEN** until a dedicated literature review; absence from this search is
not proof of novelty.

## KNOWN / CONTROL / OPEN

- **KNOWN:** vertex scaling gives the canonical tangent map `C` above.
- **KNOWN:** because the 600-cell graph is connected and contains a triangle,
  `C` is injective.  Indeed `C sigma=0` implies `sigma_u=-sigma_v` on every
  edge; an odd cycle forces one value to vanish and connectedness forces all
  values to vanish.  Hence `rank C=120` before any numerical calculation.
- **CONTROL:** the complete edge carrier and its seven minimal `2T` sectors
  were derived from literal incidence and certified independently.
- **CONTROL:** the target-free centered census and its `120:600` inertia are
  committed and byte-reproduced before this comparison.
- **OPEN:** the sign and regularity of `C^T H_p C` in every sector, variant
  and schedule.
- **OPEN:** whether `im C` is invariant under `H_p`.
- **OPEN:** whether `im C` equals the positive spectral subspace under the
  auxiliary Euclidean logarithmic-edge norm.
- **OPEN:** constraint reduction, two tensor polarizations, refinement,
  physical time, dispersion and limiting speed.

## Framing attack

Equality with the Euclidean positive eigenspace is not the invariant
definition of a conformal negative/positive mode.  An indefinite bilinear
form has many maximal definite subspaces, and its matrix eigenspaces change
under a non-orthogonal coordinate transformation.  The invariant statement
under simultaneous change of edge coordinates is instead:

```text
the restriction of the bilinear form H_p to im C is definite with the
minority sign and has the full minority dimension.
```

Therefore:

1. definiteness of `C^T H_p C`, together with the already certified full
   inertia, is the **PRIMARY** finite DeWitt-type gate;
2. invariance and equality with a matrix spectral subspace are stronger,
   coordinate-dependent **STRUCTURAL** diagnostics;
3. failure of spectral equality cannot by itself refute the primary gate;
4. a mixed or opposite-sign restricted form does refute the simple conformal
   interpretation on this fixed carrier.

There is a second limitation.  `H_p` is the Hermitian part of an averaged
off-diagonal coefficient connecting different time fibres.  It is not
already proven identical to the Lund--Regge supermetric.  Even a successful
primary gate will establish a finite action-derived DeWitt-type signature,
not a coordinate-free ghost theorem and not the existence of physical
gravitons.

## Licensed next calculation

Before reading any stored inertia eigenvector:

1. construct `C` from the literal 720 edge endpoints;
2. certify its entries, rank, `2T` equivariance and sector dimensions;
3. freeze error bounds and outcome thresholds for the restrictions
   `C_d* H_{p,d} C_d`;
4. only then inspect restricted signs and the secondary spectral angles.

No continuum harmonic, desired polarization count, dispersion relation or
speed may enter this test.
