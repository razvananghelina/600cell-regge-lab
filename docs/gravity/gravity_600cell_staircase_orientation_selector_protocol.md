# Preregistration: exact staircase-orientation selector census

Date: 2026-08-17

Prior-art gate: `7c9cd5b`.

Status: frozen before enumerating an induced cover permutation, an oriented
four-chain or the all-order orbit partition.

## 1. Frozen carrier

Rebuild the 600-cell only through `commons.cell600.build_600cell`.  Recover its
600 tetrahedra from exact adjacency incidence.  Recover the five 24-vertex
binary-tetrahedral cover cells by the same quaternion multiplication and
canonical sorting used by
`reproducible/verify_gravity_global_regge_orbits.py`.

Require before continuing:

```text
vertices = 120
edges = 720
tetrahedra = 600
cover cells = 5 x 24
every tetrahedron uses four distinct cover cells.
```

Generate all `5! = 120` total orders of the cover cells.  For each order use
the audited `build_slab` staircase formula without loading a Regge action,
metric, Hessian, canonical map or nonlinear artifact.

## 2. Complete combinatorial census

For every order require the same product-slab data:

```text
f-vector = (240, 2280, 6240, 6600, 2400)
Euler characteristic = 0
old boundary tetrahedra = 600
new boundary tetrahedra = 600
internal tetrahedral facets = 5400, each of incidence two
boundary tetrahedral facets = 1200, each of incidence one.
```

The boundary facet sets must be exactly the old and new copies of the input
600-cell.  No topological conclusion may be inferred merely from matching the
f-vector.

## 3. Frozen orientation convention

Use the stored four real vertex coordinates of the 600-cell.  For a sorted
spatial tetrahedron `(v0,v1,v2,v3)`, define its oriented-chain coefficient as
the sign of the `4 x 4` determinant whose rows are those vertex vectors.

Embed an old vertex as `(v,0)` and a new vertex as `(v,1)`.  For every sorted
staircase pentachoron define its coefficient as the sign of the `5 x 5`
determinant whose rows are these embedded vectors.  Record the smallest
absolute determinant.  A determinant with magnitude below `1e-10` is a
control failure; its sign may not be guessed.

Compute the simplicial boundary of the resulting signed four-chain exactly
over the integers.  Require:

1. every internal tetrahedral coefficient cancels to zero;
2. the old and new coefficients are opposite copies of the frozen spatial
   chain;
3. the same ordered pair of old/new signs occurs for all 120 schedules.

The sign pair may be `(-,+)` or `(+,-)` because that is a convention for the
product orientation.  It must not vary with schedule parity.

## 4. Exact symmetry action on the five cover cells

Rebuild all 14,400 `H4` vertex actions from left quaternion multiplication,
right multiplication and optional conjugation.  Retain exactly the actions
that map the five-cell cover setwise.  For each, derive its permutation of the
five cells by integer set equality.

Report:

- number of setwise actions;
- number and abstract parity census of distinct induced permutations;
- kernel size of the induced action;
- spatial orientation sign of each action, checked both from the conjugation
  construction and from its action on the signed tetrahedral chain;
- the complete orbit-size multiset on all 120 colour orders;
- for each orbit, its colour-permutation parity census and its oriented-chain
  sign pair.

Do not assume the induced group is `A5` or `S5` before enumeration.

## 5. Time reversal and flip connectivity

Let layer reversal send every vertex `v` to `v+120` and conversely.  For each
of all 120 orders, determine mechanically which schedule its reversed slab
equals.  Record whether reversal changes a symmetry orbit or colour-order
parity.

Build the graph on the 120 schedules whose edges swap one adjacent pair in
the total order.  Require each proposed edge to connect two valid slabs and
record:

- vertex and edge counts;
- connected-component sizes;
- whether every edge crosses the enumerated symmetry orbits;
- the multiset of top-simplex intersection and symmetric-difference sizes.

This is a combinatorial flip graph.  It is not evidence that the Regge action
is invariant under a flip.

## 6. Mechanical verdict

Return `ORIENTATION_SELECTS_ONE_PARITY` only if exactly one schedule orbit
passes the same frozen oriented-product and boundary-chain conditions.

Return `ORIENTATION_DOES_NOT_SELECT_PARITY` if:

1. every one of the 120 schedules passes all product-manifold and chain
   controls;
2. at least two exact symmetry orbits occur;
3. all those orbits have the same product boundary chain; and
4. the adjacent-transposition graph connects the orbits through valid
   staircase flips.

Otherwise return `OPEN_CONTROL_FAILURE`.

## 7. Interpretation boundary

A no-selector result proves only that spatial orientation, time orientation,
the relative fundamental chain and the exact `H4` cover action do not choose a
schedule parity.  It does not exclude an additional causal, variational or
matter datum.

A selector result does not validate the selected Regge dynamics; it only
removes this particular two-schedule ambiguity.

Register the verifier before its first execution.  Run only this targeted
verifier; do not run the full suite.
