# Orientation does not select a 600-cell staircase parity

Date: 2026-08-17

## Verdict

> **DERIVED COMPUTATIONAL / STRUCTURAL NEGATIVE:** spatial orientation, time
> orientation, the relative fundamental chain and the exact `H4` action on the
> fixed five-cell cover do not select between the two staircase parity classes.
> All 120 colour orders are distinct valid oriented triangulations of the same
> product slab and have the same signed boundary.

Frozen outcome:

```text
ORIENTATION_DOES_NOT_SELECT_PARITY
```

The targeted verifier passed `13/13`.  The full suite was not run.

## Provenance

- prior-art gate: `7c9cd5b`
- preregistered protocol: `5f399c2`
- verifier registered before first execution: `54f3cf4`
- verifier:
  `reproducible/verify_gravity_600cell_staircase_orientation_selector.py`
- artifact:
  `reproducible/gravity_600cell_staircase_orientation_selector.json`
- artifact SHA-256:
  `6a0bf4112baf3868a728beb45f04e9bdf1420cfbb93b26e6b8680041cb5d37f2`

No Regge action, continuum target, canonical map or nonlinear result was
loaded by the census.

## Complete finite census

The frozen five cells each contain 24 vertices, and every spatial tetrahedron
uses four distinct colours.  For each of the `5! = 120` total colour orders,
the staircase formula gives a different slab with

```text
f-vector = (240, 2280, 6240, 6600, 2400)
Euler characteristic = 0
boundary tetrahedra = 600 old + 600 new
internal tetrahedral facets = 5400, each with incidence two.
```

Thus matching counts were not used as a topology proof: exact facet incidence
and exact equality of the two boundary subcomplexes were also required.

## Oriented fundamental chains

The spatial chain was oriented by the signs of the `4 x 4` vertex-coordinate
determinants.  Each four-simplex was oriented independently by the sign of the
corresponding `5 x 5` determinant in `R^4 x R`.

The minimum absolute determinant in both calculations was

```text
0.15450849719465062,
```

far from the preregistered `1e-10` ambiguity threshold.  For every schedule,
all 5400 internal facets cancel exactly over the integers and the only
remaining chain is

```text
(old sign, new sign) = (+1,-1).
```

The overall sign is conventional.  The decisive fact is that the same
opposite pair occurs for all 120 orders, including all 60 even and all 60 odd
orders.  Permutation parity is therefore not the orientation of the
four-manifold.

## Exact `H4` cover action

The enumeration begins with all 14,400 vertex actions.  Exactly

```text
1440
```

preserve this fixed five-cell cover setwise.  They induce exactly 60 distinct
permutations of the five cells, all even, with kernel 24:

```text
induced group = A5
order orbits  = 60 even + 60 odd.
```

All 1440 setwise actions are orientation preserving, independently checked
from their quaternion-construction provenance and their action on the signed
spatial fundamental chain.  Every induced permutation transports the entire
2400-simplex slab to the predicted schedule.

This explains why the earlier complete-slab isomorphism search did not connect
the two representatives: the exact symmetry group preserves colour-order
parity.  It does not make the odd orbit invalid.

## Time reversal and flips

Layer reversal maps all 120 slabs to enumerated slabs, but it preserves the two
parity orbits:

```text
even -> even: 60
odd  -> odd:  60.
```

The adjacent-transposition graph is the connected 120-vertex, 240-edge
permutohedral graph.  Every edge crosses between the even and odd `H4` orbits.
For every edge, the two slabs have the identical census

```text
common four-simplices       1680
removed four-simplices       720
inserted four-simplices      720
symmetric difference        1440.
```

These are legitimate staircase bistellar changes.  This statement does not
assert invariance of the Regge action under them; the nonlinear calculation
has already refuted that invariance.

## Hostile framing audit

The hypothesis that orientation might rescue a unique tick is refuted on the
stated carrier.  The two nonlinear maps are not being compared across an
orientation error, a reversed time direction or an invalid triangulation.
They are two equally oriented triangulations connected inside the complete
staircase flip graph.

The census fixes one of the ten known 5-colourings of the 600-cell.  Fisk
proves the ten-colouring count in [Coloring the 600
Cell](https://arxiv.org/abs/0802.2533).  Santos' prism result explains the 120
orders and their flip graph: [Non-connected toric Hilbert
schemes](https://arxiv.org/abs/math/0204044).  The post-result search did not
locate this exact oriented five-cover schedule census.  External novelty is
therefore **OPEN**, not established by search.

## Status ledger

- **DERIVED COMPUTATIONAL:** all 120 schedules are distinct product-slab
  triangulations with the same complete oriented boundary chain.
- **DERIVED COMPUTATIONAL:** the fixed-cover `H4` action is `A5` on colours,
  with two order orbits of size 60 and kernel 24.
- **DERIVED COMPUTATIONAL:** time reversal preserves, rather than exchanges,
  those two orbits.
- **DERIVED COMPUTATIONAL:** the adjacent-transposition graph is connected and
  every one of its 240 edges crosses the two `H4` orbits.
- **STRUCTURAL NEGATIVE:** orientation and exact carrier symmetry supply no
  preferred schedule parity.
- **STRUCTURAL:** combined with the 32/32 quadratic nonlinear breaking, bare
  Regge dust on this quotient has a genuine schedule-dependent outgoing
  momentum at second order.
- **OPEN:** a selector from additional causal or matter data.
- **OPEN:** a geometry-derived improved/perfect action or suppression under a
  genuine refinement family.
- **OPEN:** behaviour on the full 720-edge phase space.

## Consequence

The current bare quotient cannot define a unique nonlinear physical tick from
the stated geometry alone.  Choosing `even` because it was enumerated first
would be fitting; choosing it because it is an even permutation would confuse
colour-order parity with manifold orientation.

The next admissible route is not another schedule search.  It is to test a
geometry-derived refinement or improved-action construction.  Without such a
result, causal speed, Planck time and nonlinear gravitational predictions from
this finite tick remain **OPEN**.
