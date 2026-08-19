# Prior-art and framing gate: static gradient versus prism-shift carrier

Date: 2026-08-19

## Question

Two independently constructed 119-dimensional spaces now occur on the same
regular 600-cell:

1. the equal-scale prism-shift carrier, represented on a tetrahedron
   `T=(v0,v1,v2,v3)` by the covector evaluations

   ```text
   a_T=(phi(v1)-phi(v0), phi(v2)-phi(v0), phi(v3)-phi(v0));
   ```

2. the static variable-connection closure kernel, represented in the local
   Cartesian frame by

   ```text
   s_T=(D_T^T)^(-1) a_T,
   D_T=[p1-p0,p2-p0,p3-p0].
   ```

Equal dimension and a common scalar parametrization do not by themselves
prove that the two embedded subspaces or their face equations coincide.  The
present target-disclosed mission asks for the literal intertwiner

```text
D_T^T s_T=a_T
```

and for equality of the two global face-matching row spaces after this
cellwise change of coordinates.

## Established mathematics

The underlying principle is standard finite-element exterior calculus.  A
continuous piecewise-affine scalar has a constant gradient on each
tetrahedron, and its tangential trace agrees across common faces.  Conversely,
on a connected tetrahedral complex with vanishing first cohomology,
tangentially matched piecewise-constant covectors are gradients modulo a
constant.  See Arnold, Falk and Winther, *Finite element exterior calculus:
from Hodge theory to numerical stability*, arXiv:0906.4325.

Therefore no novelty is claimed for the abstract isomorphism.  The only
project-specific point is whether the recently derived Poincare-translation
coordinates and spatial face transports implement exactly the same object as
the older prism-shift coordinates, without a sign, transpose or frame
convention mismatch.

## Frozen evidence

The older prism-shift calculation already establishes a 119-dimensional
potential image and exact face matching.  The newer primary and adversarial
calculations establish a 119-dimensional static closure kernel in Cartesian
translation coordinates.  Their current inputs are frozen by SHA-256 in the
subsequent protocol.

## Outcome boundary

- `SAME_CARRIER_EXACTLY` requires a cellwise invertible coordinate map,
  exact equality of potential embeddings, and equality of all global face
  row spaces.
- `SAME_ABSTRACT_SPACE_DIFFERENT_EMBEDDING` applies if the two spaces are
  isomorphic but any literal intertwining identity fails.
- `RECONCILIATION_REFUTED` applies if an upstream rank, incidence or
  potential identity fails.

A positive result is **reconciliation**, not new physics.  It closes duplicate
bookkeeping and transfers all already derived prism-shift action and canonical
elimination results to the newly named static kernel.  It supplies no time,
wave speed, graviton, mass or gravitational constant.

