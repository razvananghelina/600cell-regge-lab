# Result: global prism gluing leaves 119 longitudinal modes

Date: 2026-08-19

## Headline

The three local equal-scale prism ambiguities do not disappear when all 600
cells are glued.  Full labelled metric matching on every shared lateral face
organizes them into one exact global space:

```text
piecewise-constant local shifts with matching tangential traces
    = closed simplicial one-forms
    = gradients of vertex potentials on S3.
```

For the regular 600-cell this space has dimension

```text
120 vertices - 1 constant = 119.
```

This is not only an infinitesimal nullity.  Every vertex potential integrates
to a finite shape-matched Lorentzian polytopal slab with the same bottom,
top and vertical natural edge lengths.  Nonconstant potentials change local
four-volumes.  Therefore those lengths do not determine the static four-
geometry, even after global gluing.

The result is **DERIVED KINEMATIC, ADVERSARIALLY CORROBORATED**.  Whether the
119 modes are gauge, constraints or physical is **OPEN** until an action is
varied with respect to them.

## Provenance ledger

| stage | commit |
|---|---|
| prior-art gate | `4a9ea11` |
| frozen primary protocol | `a2c9174` |
| registered primary verifier | `165e6cb` |
| frozen primary artifact | `4d4377b` |
| adversarial protocol | `42aba25` |
| registered adversarial verifier | `5dfe65b` |
| frozen adversarial artifact | `105954e` |

Primary artifact SHA-256:

```text
1ab6654ae57c83a49dd4f427154b891c0b8ae613631773ab6733a1227b9999fa
```

Adversarial artifact SHA-256:

```text
44b0fe520c0ac46da9c0fb207bf17d72ab2611f2c1e91157122f3ad03935ea55
```

The primary verifier passes `15/15`.  The mechanically different audit
passes `14/14`.  No full suite was run.

## 1. What a shared face actually matches

Write the translation between the copies of one tetrahedron as

```text
t = N n + s,
```

where `n` is timelike normal and `s` is tangential to the spatial
tetrahedron.  On a shared triangle choose two independent edge vectors
`e1,e2`.  Its lateral triangular-prism metric sees

```text
e_i.e_j,  e_i.s,  t.t.
```

The spatial entries and `t.t=-rho` are already fixed.  Shape matching thus
equates only

```text
e1.s and e2.s.
```

It does not equate the third component normal to that triangle inside the
tetrahedron.  This distinction is load-bearing.

The rejected ambient argument had required every neighboring local shift to
be the same vector in one common `R4`.  The audit confirms that this stronger
condition has zero-dimensional intersection, but also confirms that every
face-trace map has rank two and one invisible component.  A common ambient
frame is additional structure absent from an intrinsic Regge gluing.

## 2. Exact global theorem

Let `K` be a connected nondegenerate tetrahedral complex.  On each
tetrahedron choose a constant covector.  Require its tangential restriction
to agree across every common triangular face.

Evaluating the covectors on global oriented edges gives an edge cochain
`x`.  Local constancy implies that its oriented sum around every triangle is
zero:

```text
d1 x = 0.
```

Conversely, every such edge cocycle reconstructs one constant covector on
each tetrahedron, and the shared edge values give the required face traces.
Therefore

```text
S_shift(K) is isomorphic to ker(d1),
dim S_shift(K) = V-b0+b1.
```

For a connected triangulated `S3`, `b0=1` and `b1=0`, hence

```text
dim S_shift(K)=V-1.
```

On the regular 600-cell the primary calculation formed the literal
`2400 x 1800` face-matching matrix.  Its exact rank over `F2` is `1681`.
The image of the 120 vertex potentials has exact rank `119` and lies in its
kernel.  The modular lower bound and this rational upper bound coincide, so
over the rationals and reals:

```text
ker(C)=im(B),  dim ker(C)=119.
```

The independent audit used no local matching matrix.  A spanning tree leaves
`720-120+1=601` fundamental graph cycles.  The 1,200 triangle boundaries
span all 601 over both `F3` and `F1000003`.  Five unrelated integer
potentials reconstruct exactly by tree-path integration.

## 3. Controls that could falsify the theorem

The isolated tetrahedron retains its three local modes.

The periodic `3 x 3` triangulated torus has

```text
(V,E,F)=(9,27,18),
cycle dimension=19,
triangle-boundary rank=17,
b1=2.
```

Thus its matched space has dimension

```text
V-1+b1=10,
```

not `V-1=8`.  Both primary and adversarial algorithms detect the two
topological modes.  With the triangle relations deleted, the audit also
leaves all 601 source graph cycles.  The checks therefore did not impose the
desired exactness by construction.

## 4. Finite nonuniqueness

For any vertex potential `phi` and any `rho>0`, let `G_sigma` be the positive
spatial Gram matrix of a tetrahedron and

```text
a_sigma=(phi(v1)-phi(v0),
         phi(v2)-phi(v0),
         phi(v3)-phi(v0)).
```

Define its four-dimensional cell metric by

```text
H_sigma = [ G_sigma    a_sigma ],
          [ a_sigma^T    -rho  ].
```

The Schur complement is

```text
-rho-a_sigma^T G_sigma^-1 a_sigma < 0,
```

so every cell has Lorentzian signature `(3,1)`.  Bottom and top spatial
metrics remain identical and all four struts retain squared length `-rho`.
On a common face, every mixed Gram entry is the same potential difference,
so the complete lateral face metrics match.

But the absolute determinant is

```text
|det H_sigma|
 = det(G_sigma)
   * (rho+a_sigma^T G_sigma^-1 a_sigma).
```

A nonconstant potential therefore changes the four-volume of at least one
cell without changing any declared natural edge length.

The primary one-vertex control changed 20 of 600 cell volumes, with maximum
ratio approximately `2.2197`.  The audit used squared graph distance and
`rho=7/5`; it matched all 1,200 lateral `6 x 6` interval matrices with zero
residual and changed 550 cell volumes, with maximum ratio approximately
`15.1066`.  The magnitudes are controls, not physical predictions.

## 5. Symmetry does not solve the general problem

The full declared spatial group is transitive on the 120 original vertices.
Consequently its invariant vertex potentials are constants and the coarse
600-cell has no nonzero fully invariant matched shift.  This explains why a
homogeneous coarse ansatz can consistently set shift to zero.

The canonical projected rank-edgewise carrier has ten vertex orbits of
sizes

```text
120, 600, 720, 1200, 1440,
2400, 2400, 3600, 3600, 3600.
```

It therefore has nine invariant potential gradients.  Spatial `H4`
symmetry alone does not select zero shift after refinement.

This is another reason not to identify the 119-dimensional space with the
full ADM shift.  It is a longitudinal scalar-potential sector, not three
arbitrary functions per vertex.

## 6. Consequence for the gravity route

The result closes one tempting shortcut:

> A static schedule-free cellular Regge Hessian cannot be a canonical
> function of only the bottom, top and vertical natural edge lengths.

Those data label a `V-1`-parameter family of nonisometric shape-matched
four-geometries.  Volumes already distinguish them; dihedral data need not be
single-valued functions of the declared lengths.

This does not kill Regge dynamics.  It says that the missing data must be
handled explicitly.  The next admissible object is an action

```text
S[L_bottom,L_top,rho,phi]
```

on the shape-matched family, followed by variation in `phi`.  At the
symmetric point `phi=constant`, its gradient and Hessian on the 119-
dimensional quotient decide among:

- exact null directions: candidate discrete gauge;
- nonzero constraint equations: shift multipliers or pseudo-constraints;
- nondegenerate propagation: physical internal data.

That test must precede any claim about gravitons, inertia, a limiting speed
or a physical tick.

## 7. Prior-art reconciliation

The cohomological mechanism is standard finite-element exterior calculus:
discrete de Rham complexes relate closed piecewise forms to simplicial
cohomology.  It is not a new mathematical principle:

- D. N. Arnold, R. S. Falk and R. Winther, Acta Numerica 15 (2006),
  DOI `10.1017/S0962492906210018`;
- the same authors, Bull. Amer. Math. Soc. 47 (2010), arXiv:`0906.4325`,
  DOI `10.1090/S0273-0979-10-01278-4`.

Flat-background vertex-displacement lapse/shift variables are also known in
canonical linearized Regge calculus, while curvature can break their exact
gauge symmetry:

- P. A. Hoehn, arXiv:`1411.5672`, DOI
  `10.1103/PhysRevD.91.124034`;
- B. Bahr and B. Dittrich, arXiv:`0905.1670`, DOI
  `10.1088/0264-9381/26/22/225011`.

The post-result search used the learned terms `piecewise constant closed
one-form`, `tangential trace`, `P1 gradient`, `tetrahedral prism shift
potential` and `Regge shape matching`.  It located the standard FEEC and
canonical-Regge mechanisms but no primary source for this exact
equal-scale-prism/600-cell application and its 119-mode finite family.
Search absence is not a novelty proof.  **External novelty remains OPEN.**

## 8. Status ledger

| claim | status |
|---|---|
| full face-matching kernel on the 600-cell has dimension 119 | **DERIVED EXACT** |
| it equals vertex-potential gradients modulo constants | **DERIVED EXACT** |
| the family integrates to finite Lorentzian shape-matched cells | **DERIVED KINEMATIC** |
| natural lengths determine the static cellular four-geometry | **REFUTED** |
| base `H4` symmetry admits a nonzero invariant matched shift | **REFUTED** |
| fine `H4` symmetry alone selects zero shift | **REFUTED** |
| the 119 modes are gauge | **OPEN** |
| the 119 modes are the full ADM shift | **NOT ESTABLISHED / DIMENSIONALLY DIFFERENT** |
| a shift-extended cellular action is canonical | **OPEN** |
| external novelty | **OPEN** |

Only the two mission verifiers and static guards are run.  The full suite is
not run by explicit user instruction.

