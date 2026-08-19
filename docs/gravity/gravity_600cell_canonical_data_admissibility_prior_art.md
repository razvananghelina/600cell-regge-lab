# Prior-art and framing gate: schedule-free canonical-data admissibility

Date: 2026-08-19

## Exact question

Fix the regular 600-cell lower slice and one of the already certified
nonstatic homothetic backgrounds

```text
q_i=lambda p_i+tau n,  (lambda,tau)=(2,5) or (3,11).
```

Each of the 600 spacetime cells is required to remain one flat Lorentzian
tetrahedral frustum.  Across every one of the 1200 lateral triangular-prism
faces, allow exactly the already derived one-dimensional variation of the
Poincare transition fixing the lower triangle pointwise.  The prospective
canonical data, with the lower slice fixed, are

```text
720 independent upper spatial squared-edge variations
+120 independent corresponding-strut squared-length variations
=840 data coordinates.
```

The question is whether all 840 data directions admit compatible local
frustum deformations, or only a proper subspace does.  Equivalently, after
the six fixed-length flex coordinates of every cell are included and the
face-transition scalar is eliminated, what is the projection to the 840
global data coordinates of the complete linearized gluing kernel?

This is a kinematic tangent-space question.  It does not assume or evaluate
an action, Hessian, equation of motion, lapse selection, graviton, continuum
limit or speed.

## Why this gate is necessary

The existing complete variable-face matrix has shape `6000 x 3600`.  At the
two nonstatic representatives it has exact rational column rank 3600,
certified by nonzero minors modulo two primes.  That proves uniqueness of a
compatible flex when all natural lengths are held fixed.  It does **not**
prove existence for an arbitrary inhomogeneous variation of those lengths.

Injectivity of the homogeneous operator is not surjectivity of the forced
problem.  Promoting it to a full canonical carrier without adding the forcing
columns would be a linear-algebra error.  In particular, a simplicial
staircase always supplies an action for its own artificial diagonal data,
but that fact does not prove that a single flat schedule-free frustum exists
for arbitrary boundary and strut data.

Therefore the proposed refinement/pseudo-constraint Hessian is not yet
authorized.  Its domain must first be measured rather than assumed.

## What is already established in the repository

- **DERIVED EXACT:** one isolated frustum has a rank-ten natural-length map
  in sixteen upper-vertex coordinates and hence a six-dimensional fixed-data
  kernel.
- **DERIVED EXACT:** allowing the unique lower-face stabilizer in a glued
  pair changes the compatible dimension from six to seven and supplies one
  relative connection-coupled direction.
- **DERIVED EXACT:** after eliminating that face scalar, every face contributes
  a rank-five block on the two local six-dimensional flex spaces.
- **DERIVED EXACT:** on the complete 600-cell the nonstatic fixed-data matrix
  has rank 3600/3600, whereas the static matrix has the already reconciled
  119-dimensional P1-gradient kernel.
- **NOT COMPUTED:** the columns induced by independent upper-edge and strut
  variations, their augmented rank, and the dimension of their admissible
  projection.

The new calculation must not rediagonalize or reinterpret the static 119
modes.

## Primary literature boundary

Connection and shape-matching formulations already explain why local frame
or connection data are larger than ordinary length-Regge data until
metricity/gluing constraints are imposed:

- Dittrich and Ryan, *Phase space descriptions for simplicial 4d
  geometries*, arXiv:[0807.2806](https://arxiv.org/abs/0807.2806), separate
  connection/area data from the length-Regge sector and impose gluing and
  simplicity restrictions.
- Anza and Speziale, *A note on the secondary simplicity constraints in loop
  quantum gravity*, arXiv:[1409.0836](https://arxiv.org/abs/1409.0836), show
  that Lorentzian secondary constraints recover Levi-Civita and shape-
  matching information only on the constrained sector.
- Khatsymovsky, *Affine connection form of Regge calculus*,
  arXiv:[1509.04974](https://arxiv.org/abs/1509.04974), treats affine
  transition matrices as independent variables and recovers ordinary Regge
  calculus only after their equations are imposed.
- Bahr and Dittrich, *[Broken Gauge Symmetries and Constraints in Regge
  Calculus](https://arxiv.org/abs/0905.1670)*, establish that curved discrete
  backgrounds generically replace exact gauge constraints by background-
  dependent pseudo-constraints.
- Hoehn, *[Canonical linearized Regge Calculus: counting lattice gravitons
  with Pachner moves](https://arxiv.org/abs/1411.5672)*, identifies exact
  vertex-displacement constraints and curvature-carrying lattice gravitons
  on flat linearized simplicial backgrounds.

These sources make the distinction between natural data and constrained
connection data standard.  They do not calculate the present 600-cell
augmented compatibility rank or prove that all 840 directions are admitted.

## Classification before calculation

- **KNOWN:** length, connection and shape-matching variables need not have
  the same unconstrained tangent space.
- **CONTROL:** the old fixed-data columns must retain rank 3600 at both
  nonstatic representatives; the local natural-length Jacobian must retain
  rank ten.
- **OPEN:** the exact augmented modular rank and admissible-data dimension.
- **OPEN:** whether the admissible subspace is all 840 directions, only the
  homogeneous scale/lapse plane, or an intermediate carrier.
- **NOT CLAIMED:** external novelty.  A targeted search found no source for
  this exact finite matrix, but search absence is not a novelty proof.

## Consequences fixed in advance

1. If all 840 directions are admitted with an exact certificate, the
   schedule-free cellular geometry passes the tangent-domain gate.  This
   authorizes a separate action/Hessian protocol; it does not itself select a
   tick or produce a wave speed.
2. If only a proper subspace is admitted, its dimension and explicit content
   must be classified before any action is restricted to it.  A large
   schedule-dependent Hessian remains unauthorized.
3. If only modular deficits are obtained without rational witnesses, the
   dimension remains **OPEN**; agreement at two primes is not promoted to a
   rational nullity theorem.
4. A zero admissible dimension is impossible if the implementation is
   correct, because the two-parameter homothetic family supplies scale and
   strut controls.  Failure to recover those controls is an implementation
   failure, not a negative physical result.

The calculation is therefore a domain-of-definition test for the next
gravity step, not evidence for gravity, time, `c` or `G` by itself.
