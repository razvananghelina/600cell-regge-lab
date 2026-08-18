# Prior-art gate: oriented orbits of all ten 600-cell colourings

Date: 2026-08-17

Status: written before computing the chromatic degree of the other nine
five-colourings or their orientation-preserving symmetry orbits.

## Exact object, carrier and hypotheses

Let `K` be the determinant-oriented boundary complex of the regular 600-cell
on the 120 quaternionic vertices reconstructed by `commons/cell600.py`.

A colouring is an *unordered* partition of the 120 vertices into five
independent 24-sets.  It is not generated from one favoured left- or
right-coset construction: the census must first prove that all maximum
independent 24-sets and all exact five-set covers have been included.

For each colouring `C` and every total order of its five cells, use the order
as the labels `0,...,4`.  The resulting simplicial map

```text
f_(C,order) : K -> boundary(Delta^4)
```

has the integral chromatic degree defined in the previous fixed-cover audit.
Call an ordered colouring *orientation-compatible* when this degree is
positive relative to the frozen determinant orientation of `K` and the
standard orientation `(0,1,2,3,4)` of the target simplex.  This sign convention
is target-blind: no Regge action or schedule output is loaded.

The group acting on compatible ordered colourings is the index-two,
orientation-preserving subgroup `H4+` of order 7200, realized by

```text
q -> l q r^-1.
```

The conjugating maps `q -> l conjugate(q) r^-1` reverse the ambient
four-dimensional orientation and are retained only as an improper-symmetry
control.  They are not allowed to identify two configurations of a fixed
oriented carrier.

## Primary prior art

- Steve Fisk proves that the 600-cell has exactly ten five-colourings and
  describes the left- and right-coset constructions: [Coloring the 600
  Cell](https://arxiv.org/abs/0802.2533).
- The ten decompositions into five disjoint 24-cells go back to Schoute and
  are recorded by Coxeter.  A modern explicit treatment of the 25 embedded
  24-cells through left and right binary-tetrahedral cosets is Theorem 4.18
  of [Geometry and combinatorics of the 600-cell and the
  120-cell](https://repository.tudelft.nl/file/File_fd61d63c-22a4-46f7-8992-cf7b511df139).
- The combinatorial degree of a colouring map to a simplex boundary and its
  sign change under an odd target relabelling are standard; see Definition
  1.11 of [Toric residue and combinatorial
  degree](https://arxiv.org/abs/math/0309409).

The located primary sources establish the ten colourings and the degree
formalism.  They do not report the integer degree of all ten 600-cell
colourings or the orbit census of positive-degree ordered colourings under
`H4+`.

## KNOWN

- The repository already has a solver-certified, target-blind census of 25
  maximum independent 24-sets and exactly ten exact covers.  The full order-
  14400 `H4` group is transitive on the ten unordered covers.
- The full `H4` action splits all `10*120 = 1200` ordered covers into two
  600-element orbits with pointwise stabilizer 24.
- For one fixed cover, the chromatic degree is nonzero with magnitude 72 and
  its 120 orders split into 60 of each sign.
- Proper rotations preserve the source orientation; improper symmetries
  reverse it.

## Disclosed prediction before the calculation

Orbit-stabilizer counting suggests, but does not yet prove, the following:

```text
7200 / 24 = 300,
```

so the 600 positive-degree ordered colourings may split into two proper-
rotation orbits of size 300.  Equivalently, the ten unordered covers may split
into two chiral proper-rotation orbits of five, exchanged only by an improper
symmetry.  This expectation is disclosed now and cannot later be presented as
a blind discovery.

The calculation is allowed to refute it.  In particular, stabilizers might
contain orientation-reversing elements, degrees might differ between covers,
or the compatible set might have a different orbit structure.

## CONTROL

- Reconstruct the 25 candidate cells as the full `H4` orbit of `2T`.
- Prove `alpha=24` and exclude any further independent 24-set by exact CP-SAT,
  then enumerate the exact covers.  This makes the ten-cover input
  exhaustive rather than assumed.
- Compute each cover degree independently from all five target facets and
  direct signed preimage counts.
- Enumerate all 120 orders per cover and require exact alternation with order
  parity.
- Construct all 7200 proper and all 7200 improper vertex permutations,
  verify their determinant/orientation class independently, and compute the
  complete orbit/stabilizer census.
- Load no Regge action, nonlinear schedule result, preferred parity,
  continuum target or physical scale.

## OPEN and decision boundary

- **OPEN:** whether all ten degrees have magnitude 72.
- **OPEN:** whether compatible ordered colourings form one proper-rotation
  orbit, two chiral orbits, or more.
- **OPEN:** whether any physical axiom requires chromatic compatibility at
  all.

Mechanical outcomes:

- `ONE_ORIENTED_CANONICAL_CLASS`: all compatible ordered colourings form one
  `H4+` orbit.
- `CHIRAL_COVER_AMBIGUITY`: all controls pass but more than one `H4+` orbit
  remains, with improper symmetries relating at least two of them.
- `OPEN_CONTROL_FAILURE`: the cover census, degree controls or symmetry
  controls fail.

Even the first outcome would establish only a **STRUCTURAL** canonical class.
It would not turn chromatic orientation into a physical time axiom.  The
external novelty of any finite orbit result remains **OPEN** pending a
dedicated literature review.
