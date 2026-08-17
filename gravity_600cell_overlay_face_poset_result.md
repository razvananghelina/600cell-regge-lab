# Universal overlay face poset and global order complex

Date: 2026-08-17

Only the targeted exact face-poset verifier was run.  No metric, gravity
action or full suite was run.

## 1. Provenance

- prior-art gate: `4a85d25`
- preregistered protocol: `e8d995c`
- verifier registered before enumeration: `a139cfa`
- verifier:
  `reproducible/verify_gravity_600cell_overlay_face_poset.py`
- frozen overlay source SHA-256:
  `0dd03eed878f599463a44160484c74ddeaa0511fc70c8b2e77bc05a2f36dd3dc`
- result artifact:
  `reproducible/gravity_600cell_overlay_face_poset.json`
- result artifact SHA-256:
  `439a3d067d50415f0a47c79091ec746c12dd7975b2246b6143f3f7a70847ce13`

The verifier passed `13/13` controls using exact rational intersections,
affine ranks and face incidences.

## 2. Exact local polyhedral complex

Of the `binomial(20,4)=4,845` possible four-wall subsets, 3,282 have
independent affine equations.  After exact prism containment and rational
deduplication, the arrangement has 33 vertices.

The complete local polyhedral f-vector is

```text
(f0,f1,f2,f3,f4) = (33,206,468,442,148),
Euler characteristic = 1.
```

Its boundary has

```text
(22,100,152,74),
Euler characteristic = 0.
```

Each time end is the unsplit tetrahedral face lattice `(4,6,4,1)`.  Every
interior three-face belongs to exactly two four-cells and every boundary
three-face to exactly one.  Every positive-dimensional cell independently
passes the Euler test for the sphere bounding a convex polytope.

These facts certify a regular polyhedral four-ball.  They are **DERIVED
COMBINATORIAL**.

## 3. Exact symmetry

All 48 elements of `S4 x C2` preserve every rational vertex, face, dimension,
inclusion and boundary label.  The 1,297 faces form 80 symmetry orbits,
distributed by dimension as

```text
dimension 0:  6 orbits
dimension 1: 15 orbits
dimension 2: 22 orbits
dimension 3: 23 orbits
dimension 4: 14 orbits.
```

Thus the carrier can be represented locally by a modest number of orbit
types, even though its global expansion is large.

## 4. Canonical local barycentric subdivision

The order complex has one vertex per polyhedral face and one simplex per
strict face chain.  Its exact local f-vector is

```text
(1297,14436,44004,50832,19968),
Euler characteristic = 1.
```

Its boundary vector is

```text
(348,2172,3648,1824),
Euler characteristic = 0.
```

Each time end gives the barycentric tetrahedron vector `(15,50,60,24)`.
Because the order-complex construction is functorial, adjacent prism
triangulations agree without a new diagonal choice.

## 5. Exact global 600-cell assembly

Classifying every face chain by the minimal spatial simplex supporting its
maximal face and assembling with the certified spatial f-vector
`(120,720,1200,600)` gives the global polyhedral vector

```text
(9960,86040,230880,243600,88800),
Euler characteristic = 0,
```

and the global barycentric simplicial vector

```text
(659280,7977360,25303680,29966400,11980800),
Euler characteristic = 0.
```

Both time boundaries reproduce exactly

```text
(2640,17040,28800,14400),
```

the already certified barycentric 600-cell boundary carrier.  This is a
strong independent gluing control.

## 6. Size verdict

The carrier is exact but expensive:

- one ordinary staircase: `600*4 = 2,400` four-simplices;
- naive product barycentric carrier: 115,200 four-simplices;
- universal-overlay order complex: 11,980,800 four-simplices.

Thus the universal carrier is 4,992 times larger than one staircase and 104
times larger than the naive barycentric product.  It also contains 7,977,360
edges, which would be the unconstrained metric-variable count in a completely
generic edge-length Regge model.

Therefore:

- **DERIVED POSITIVE:** the canonical common simplicial carrier exists and
  glues globally;
- **DERIVED PRACTICAL NEGATIVE:** a generic direct Regge evolution on all its
  edge variables is outside the present computational route;
- **STRUCTURAL POSITIVE:** a homogeneous calculation can be performed on the
  19,968 local four-simplices and then assembled by symmetry, so the next
  Lorentzian gate remains feasible.

## 7. Post-result prior art

The June 2026 paper
[*Chamber geometry and specification numbers of Boolean threshold
functions*](https://arxiv.org/abs/2606.29477) identifies threshold functions
with chambers of the threshold arrangement and essential Boolean points with
chamber facets.  Hence the chamber/facet interpretation is **KNOWN**.  The
specific bounded-sector f-vector and its 600-cell order-complex assembly were
not located there; external novelty of those counts remains **OPEN**.

## 8. Status ledger

| Claim | Status |
|---|---|
| Local polyhedral f-vector `(33,206,468,442,148)` | **DERIVED EXACT** |
| Local complex is a four-ball with correct boundary incidences | **DERIVED EXACT** |
| Full `S4 x C2` face-poset symmetry | **DERIVED EXACT** |
| Global polyhedral f-vector | **DERIVED EXACT** |
| Global order-complex f-vector and 11,980,800 top simplices | **DERIVED EXACT** |
| Both time boundaries equal the certified 2640 carrier | **DERIVED EXACT** |
| Generic full-edge evolution is practical now | **DERIVED PRACTICAL NEGATIVE** |
| Homogeneous Lorentzian local realization | **OPEN** |
| Every refined simplex is nondegenerate on an evolving slab | **OPEN** |
| Regge--dust action on the refined carrier | **OPEN** |
| Continuum or schedule-independent dynamics | **OPEN** |

## 9. Next falsifier

The next computation must not materialize 12 million simplices.  It should
enumerate the 19,968 local maximal face flags and assign their barycentric
vertices physical coordinates from a completely stated homothetic slab map.
It must then compute exact determinant polynomials in the lower scale, upper
scale and lapse and decide whether any simplex becomes degenerate on the
physical positive branch.

At equal lower and upper scales the product map is affine and nondegeneracy is
mostly a control.  The load-bearing case is unequal scales, where the
homothetic interpolation and the choice of physical face barycentres must be
justified before calculation.  If that choice is not canonical, the refined
Regge route stops here rather than hiding a new scheduling freedom.

