# Prior-art gate: direct-precision action-York residual

Date: 2026-08-19

## Question disclosed before the calculation

The preceding target-disclosed calculation left the exact identity

```text
generalized negative shape carrier = action-weighted tangential image
```

OPEN.  In all sixteen selected cells its binary midpoints gave

```text
||P_L-P_-||_2 about 4.07209629e-4,
||L* A T||_2 about 2.96331560e-5,
```

but the entrywise Flint-ball propagation through the two Legendre inversions
was too broad to classify either quantity as resolved nonzero.  These target
values are known before the present protocol; this is not a blind search.

The direct-precision question is whether the second value survives a fresh
construction from the local Regge Hessians before the interval-wrapping step.

## Literature search

Searches performed after the OPEN result included:

```text
Regge calculus curved background vertex displacement symmetry broken pseudo constraints
linearized Regge calculus curved background gauge modes Hessian vertex displacement
discrete gravity pseudo constraints broken diffeomorphism symmetry curved Regge
```

The closest primary sources are:

1. Bahr and Dittrich, [(Broken) Gauge Symmetries and Constraints in Regge
   Calculus](https://arxiv.org/abs/0905.1670).  Curved Regge solutions do not
   generically retain exact discrete gauge symmetries; canonical constraints
   are replaced by pseudo-constraints.
2. Dittrich and Hoehn, [From covariant to canonical formulations of discrete
   gravity](https://arxiv.org/abs/0912.1817).  Exact vertex-displacement
   constraints occur in the linearized flat-background theory, while higher
   orders make them background-dependent pseudo-constraints.
3. Hoehn, [Canonical linearized Regge Calculus: counting lattice gravitons
   with Pachner moves](https://arxiv.org/abs/1411.5672).  The flat linearized
   theory separates vertex-displacement generators from gauge-invariant
   lattice gravitons.

No source found computes the action-weighted York carrier of the fixed
600-cell dust slab, its two one-dimensional binary-tetrahedral sectors, or the
specific residual above.  External novelty is OPEN.

## Framing attack

A stable nonzero midpoint at one approximate background is not automatically
a theorem about an exact continuum gauge identity.  Three issues must be kept
separate:

1. **Interval wrapping.**  A broad enclosure produced after repeated matrix
   inversions can fail to resolve a genuinely nonzero residual.  Repeating the
   same enclosure route cannot settle this.
2. **Background accuracy.**  A computation at an uncertified approximate
   solution could split an identity that holds at the exact solution.  The
   frozen first and second tick artifacts therefore enter with their complete
   state strings, junction bounds and residuals; the direct calculation must
   report those controls.
3. **Discrete versus continuum symmetry.**  Resolving a nonzero residual
   refutes the exact identity on this fixed curved discrete slab.  It does not
   refute continuum diffeomorphism symmetry.  Conversely, resolving zero on
   this slab would be STRUCTURAL, not a continuum gauge theorem.

The learned quantity to test is basis-independent generalized invariance.  If
`L` spans the geometry-selected longitudinal image, it must satisfy

```text
A L subset B L
```

to equal a generalized eigenspace of `A e = lambda B e`.  This can be tested
directly by the residual of `A L` from `im(B L)`, in addition to the earlier
cross block `L* A T`.  Either resolved nonzero residual falsifies exact
equality without relying on eigenvalue labels or dimension matching.
