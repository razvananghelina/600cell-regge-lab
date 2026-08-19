# Prior-art gate: canonical elimination of the relative prism shifts

Date: 2026-08-19

Status: **target-disclosed composition audit before its protocol and
verifier**.

## 1. Inputs already derived

Two independent repository results now meet on the same 120 vertical-edge
coordinates.

First, the exact unequal-scale frustum theorem gives

```text
(q-1)*s=0
```

for common struts, and for nonuniform strut squares gives an invertible local
reconstruction

```text
s=E^-1*c/[2*(q-1)].
```

On the 600-cell, relative vertex-strut variations have dimension `119`; the
one complementary direction is the collective lapse.

Second, the complete non-static dust-slab calculation partitions the
canonical pre-Legendre Jacobian as

```text
J=[A B],
  [C D]
```

with 120 pole columns/equations and 1,440 strong columns/equations.  It has
already certified

```text
A invertible,
S_pole=D-C*A^-1*B invertible, rank 120/120
```

in both staircase parities and all seven binary-tetrahedral sectors.

The new question is not another spectrum search.  It is whether these two
facts compose to classify the 119 relative shift/strut modes as free
propagating data or as variables fixed by the canonical equations.

## 2. Complete hypotheses

The proposed conclusion is conditional on all of the following:

1. the accepted first non-static fixed-mass dust tick;
2. either of the two already derived order-24 staircase parities;
3. the complete `2,280` logarithmic signed-squared edge carrier;
4. fixed old-boundary geometry and fixed incoming canonical momentum;
5. the `1,560 x 1,560` pre-Legendre linearization with 840 internal and 720
   new-boundary unknowns;
6. the 120 vertical pole magnitudes as the weak coordinates `z`;
7. all other internal diagonals and all new-boundary variables in the strong
   coordinate `x`;
8. the exact 119-dimensional relative-pole subspace, with the collective
   pole direction excluded;
9. only a local linear statement at the accepted background.

This is not a statement about the equal-scale static action alone, arbitrary
nonhomogeneous boundaries, nonlinear global uniqueness, a continuum gauge
quotient, tensor gravitons or a refinement limit.

## 3. Framing attack: the geometric embedding changes diagonals

The exact polytopal reconstruction does not vary a pole in isolation.  A
relative strut variation also changes the cross-diagonal lengths needed to
realize the tangential shift.  Therefore it would be invalid merely to point
at the pole block `D` and call the mode stiff.

However, those diagonal variations form a graph inside the strong variables:

```text
x=G*z+y.
```

Under this change the linear equations become

```text
[A, A*G+B],
[C, C*G+D].
```

Eliminating `y` gives

```text
(C*G+D)-C*A^-1*(A*G+B)=D-C*A^-1*B=S_pole.
```

Thus the effective pole equation is exactly independent of which legitimate
strong-coordinate graph represents the geometric frustum.  This argument
requires `A` to be invertible and the geometric diagonal changes to belong to
the declared strong carrier; both conditions must be mechanically checked.

If `S_pole` and the relative embedding are injective, the homogeneous
canonical equations imply

```text
S_pole*z=0  =>  z=0.
```

That would refute a free propagated `phi` field.  It would not prove that the
variables have no sourced boundary effect, nor identify the canonical lift
with the naive geometric graph.

## 4. Prior art

The covariant-to-canonical Regge framework and its background-dependent
pseudo-constraints are established in:

- B. Dittrich and P. A. Hoehn, *From covariant to canonical formulations of
  discrete gravity*, [arXiv:0912.1817](https://arxiv.org/abs/0912.1817);
- B. Bahr and B. Dittrich, *(Broken) Gauge Symmetries and Constraints in
  Regge Calculus*, [arXiv:0905.1670](https://arxiv.org/abs/0905.1670);
- P. A. Hoehn, *Canonical linearized Regge Calculus: counting lattice
  gravitons with Pachner moves*,
  [arXiv:1411.5672](https://arxiv.org/abs/1411.5672).

These sources explain why a curved discrete background can replace exact
lapse/shift constraints by nonzero pseudo-constraints.  The Schur-complement
identity itself is elementary linear algebra and carries no novelty claim.

The focused search used `Regge bulk Hessian Schur complement`, `integrate out
internal edges`, `boundary effective action` and `pseudo-constraint Hessian`.
It located the general canonical mechanism above, but no source composing the
present exact 600-cell frustum branch with the certified 120-pole Schur
operator.  Search absence is not proof; external novelty remains **OPEN**.

## 5. What can and cannot be concluded

If the composition passes:

- **DERIVED:** the 119 relative pole/shift coordinates are not free solutions
  of the fixed-data linearized canonical equations;
- **DERIVED:** the conclusion is invariant under the strong diagonal graph
  chosen to realize the polytopal shift;
- **STRUCTURAL:** the sector is auxiliary/pseudo-constraint-like on this
  curved dust background;
- **OPEN:** whether a sourced response transmits into physical boundary
  tensor modes after constraint reduction.

If any Schur sector is zero/open, or if the relative map fails rank 119, the
elimination conclusion is **OPEN**.  No desired dispersion or speed may be
used to repair it.
