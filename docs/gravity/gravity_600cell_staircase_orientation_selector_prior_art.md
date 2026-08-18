# Prior-art gate: can orientation select a 600-cell staircase parity?

Date: 2026-08-17

Status: written before enumerating all phase orders, computing the induced
five-cell permutation group or evaluating an oriented fundamental chain.

## Exact object and hypotheses

Let `K` be the fixed triangulated 600-cell boundary with 120 vertices, 720
edges and 600 tetrahedra.  Retain the already derived partition of its
vertices into five independent binary-tetrahedral cover cells of size 24.
Every spatial tetrahedron contains four distinct cover colours.

For a total order `sigma` of the five colours, triangulate each product prism

```text
tetrahedron x [0,1]
```

by the same staircase rule used in
`reproducible/verify_gravity_global_regge_orbits.py`.  Denote the resulting
2400-pentachoron slab by `T_sigma`.  The boundary vertices and the old-to-new
vertical pairing are fixed pointwise; only the internal staircase diagonals
change.

The question is not whether the even and odd Regge actions differ.  That is
already derived.  The question is whether the orientation of `K`, the
orientation of `[0,1]`, the induced relative fundamental class of
`K x [0,1]`, or an exact `H4` action makes only one parity of `sigma`
admissible.  If so, comparing the two parities as two physical schedules was
too broad.  If both parities triangulate the same oriented product with the
same oriented boundary, orientation cannot select the dynamics.

No Regge action, Hessian, nonlinear output, continuum value or desired parity
will be loaded in the enumeration.

## Primary prior art

Santos proves that triangulations of a simplex prism `Delta_d x I` are in
bijection with linear orders of its vertices and that adjacent transpositions
give bistellar flips.  He also identifies compatible refinements of
`T x I` with locally acyclic orientations of the one-skeleton:
[Non-connected toric Hilbert schemes](https://arxiv.org/abs/math/0204044),
Propositions 1.2--1.3.

The prism/permutation and flip correspondence is also treated in
[Triangulations of prisms and preprojective algebras of type
A](https://arxiv.org/abs/2208.12957).  General geometric bistellar flips are
reviewed in [Geometric bistellar flips: the setting, the context and a
construction](https://arxiv.org/abs/math/0601746).

For the gravitational interpretation, standard Regge dynamics generally
remembers triangulation at nonlinear order; exact discrete symmetry is not
automatic on curved solutions:

- [(Broken) Gauge Symmetries and Constraints in Regge
  Calculus](https://arxiv.org/abs/0905.1670);
- [From covariant to canonical formulations of discrete
  gravity](https://arxiv.org/abs/0912.1817);
- [Improved and Perfect Actions in Discrete
  Gravity](https://arxiv.org/abs/0907.4323).

## KNOWN / CONTROL / OPEN

### KNOWN

- A total colour order supplies a globally compatible locally acyclic
  orientation and hence a staircase triangulation of the product slab.
- Adjacent colour transposition is a legitimate bistellar change, not by
  itself a reversal of time or of the manifold orientation.
- The already enumerated even and odd representatives are not related by a
  complete `H4` slab isomorphism, while their linear canonical boundary maps
  coincide and their nonlinear maps differ quadratically.

### CONTROL

- Rebuild the 600-cell, cover and all 120 colour orders from exact incidence.
- Require every order to have 2400 distinct four-simplices, the fixed two
  600-cell boundary complexes, manifold codimension-one incidence and the
  product Euler characteristic.
- Orient the spatial tetrahedra coherently, then orient every pentachoron by
  the product embedding in `R^4 x R`; require cancellation on every internal
  tetrahedral facet and the same signed boundary chain
  `K_new - K_old` for every accepted order.
- Enumerate the setwise `H4` stabilizer of the five-cell cover and its exact
  induced permutation group.  Compute its orbits on all 120 orders without
  assuming that they are the two parity classes.
- Build the adjacent-transposition graph on all orders and verify which
  symmetry orbits it connects.

### OPEN

- Whether all 120 global orders pass the oriented-product controls.
- Whether the induced cover action is exactly `A5`, `S5` or a smaller group.
- Whether the 120 orders form one or more exact `H4` orbits.
- Whether either orbit is excluded by the fixed spatial and temporal
  orientation.

## Framing attack and decision rule

It is not enough to observe that the two representatives have opposite
permutation parity.  Permutation parity is a label of the chosen colour order,
not automatically the orientation of the triangulated four-manifold.

Call `ORIENTATION_SELECTS_ONE_PARITY` only if exactly one symmetry orbit gives
a coherent relative fundamental chain with the frozen boundary and time
orientation.  Call `ORIENTATION_DOES_NOT_SELECT_PARITY` if both schedule
orbits pass the same oriented-product chain and are linked by legitimate
staircase flips.  Otherwise return `OPEN_CONTROL_FAILURE`.

The second outcome would sharpen the nonlinear no-go: both choices are valid
oriented triangulations, so bare Regge dynamics needs an additional selector,
a perfect action or a refinement limit.  It would not prove that no other
geometric datum can select a schedule.

External novelty of the exact 600-cell census remains **OPEN**.
