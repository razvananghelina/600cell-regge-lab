# Adversarial protocol: ordered-colour ambiguity of the balanced slab

Date: 2026-08-19

Primary protocol commit: `f39a5cc`.
Primary verifier commit: `115fef4`.
Primary artifact commit: `fc18efd`.

Primary artifact SHA-256:

```text
0a9e9e796cd671c82f2e428bfa21ba63ccb07fe76867e4553979c3c54b22a0d5
```

## Result under attack

The primary route found that the rank-edgewise complex has a proper derived
four-colouring and a conforming staircase product, but that the four colour
classes admit `4!=24` distinct linear orders.  All 24 are spatially
`H4`-invariant; time reversal pairs them into twelve orbits.  Under the frozen
selection rule this fails canonical selection.

The audit must attack the decisive count and the interpretation, not rerun
the full carrier code.

## Independent mechanism

The audit will not construct the complete 460,800-pentachoron slab and will
not propagate colours through the 115,200-tetrahedron dual graph.

Instead it will:

1. reconstruct one fine ranked chamber from the Edelsbrunner--Grayson
   weak-composition colour schemes rather than the primary eight-child list;
2. build its graph with NetworkX and enumerate proper four-colourings by an
   independent exact backtracker;
3. encode a prism schedule only by the six square-face diagonals
   `edge x I`.  Each linear colour order orients all six pairs and hence gives
   a tournament on four colours;
4. prove computationally that the admissible tournaments are exactly the 24
   transitive tournaments.  Because the six boundary-square diagonals can be
   recovered from the slab, different tournaments give different global slab
   complexes on any carrier containing a tetrahedron;
5. reverse every tournament and count the resulting orbits;
6. test whether spatial rank-preserving actions can identify any two orders.

This is mechanically different from comparing full pentachoron sets.

## Interpretation attack

Three input structures will be separated explicitly:

```text
D0: uncoloured carrier;
D1: carrier plus an unordered partition into four colour classes;
D2: carrier plus a map into the linearly ordered set 0<1<2<3.
```

`D2` selects one staircase by definition.  The audit must not misreport that
tautology as geometric selection.  The question is whether the already
derived carrier supplies `D2`, or whether declaring the linear order on the
modulo-four colour classes is extra structure.  The primary protocol froze
the latter reading; the audit cannot change that criterion retroactively.

The post-result literature check confirms that Joswig--Witte take linear
vertex orders as input and explicitly warn that different orders can give
non-isomorphic product triangulations.  Their theorem establishes conformity,
not selection.  No searched source proves that the modular colouring of an
edgewise subdivision carries a geometrically forced linear order for a
temporal product.  Search absence is not a no-go theorem, so this exact
interpretive question remains separately labelled.

## Controls

- **positive control:** a segment times an interval has two staircase
  diagonals, exchanged by time reversal;
- **negative control:** a directed three-cycle inside a four-colour
  tournament is not induced by any total order and must be rejected;
- **orientation attack:** even if a spatial orientation selected permutation
  parity, the 24 orders could shrink only to two classes of 12, not to one.

## Frozen verdict rule

- If the independent diagonal/tournament encoding does not give 24 distinct
  schedules, the primary result is **OPEN** and the discrepancy is the result.
- If it gives 24, the combinatorial ambiguity is adversarially corroborated.
- If an input structure already certified before this mission supplies an
  ordered colour map `D2` without a new convention, that fact may challenge
  the framing but does not retroactively pass the preregistered gate.
- Otherwise the slab exists but is **STRUCTURAL**, and the local-lapse mission
  stops before choosing a schedule.

The audit will be registered and only the two mission verifiers plus static
guards will be run.  No full suite is allowed.

