# Result: balanced slabs exist, but temporal order is not selected

Date: 2026-08-19

## Headline

Both certified projected spatial carriers admit an exact global proper
four-colouring and a conforming staircase product with an interval.  On the
fine carrier the standard weighted-rank residue is

```text
c(v,w)=r(v)+r(w) mod 4.
```

The existence result is clean.  The selection result is negative: the four
colour classes admit `4!=24` distinct staircase orders, every one is
compatible with the full tested spatial `H4` action, and time reversal only
pairs them into twelve size-two orbits.  An independent tournament audit
recovers the same count without constructing the full slab.

Under the preregistered boundary, the temporal carrier is **STRUCTURAL**, not
selected.  No local lapse or Hessian was computed on one chosen schedule.

## Provenance ledger

| stage | commit |
|---|---|
| prior-art gate | `dabc098` |
| primary protocol | `f39a5cc` |
| registered primary verifier | `115fef4` |
| frozen primary artifact | `fc18efd` |
| adversarial protocol | `8e4c2fb` |
| registered tournament audit | `50feab0` |
| frozen adversarial artifact | `5645d2b` |

The primary artifact was reproduced byte-for-byte with SHA-256

```text
0a9e9e796cd671c82f2e428bfa21ba63ccb07fe76867e4553979c3c54b22a0d5.
```

The adversarial artifact was reproduced byte-for-byte with SHA-256

```text
dd1043a8cb712adb4f0717f95024b9ce62132501198938bb997e7ab3dad8bf65.
```

## 1. Spatial reconstruction and colouring

The primary verifier rebuilt the 600-cell from adjacency, formed all
barycentric flags and then used the direct eight-child rank split.  It did
not import a previous slab or a Regge action.

The two spatial f-vectors were reproduced exactly:

```text
P(sd K_600)             (2640, 17040,  28800,  14400)
P(Esd_2(sd K_600))     (19680,134880, 230400, 115200).
```

Every triangle has incidence two and both tetrahedron dual graphs are
connected.  The base face-rank colouring has class sizes

```text
(120,720,1200,600),
```

and the fine weighted-rank colouring has class sizes

```text
(4920,3840,4920,6000).
```

Every spatial edge has differently coloured endpoints and every tetrahedron
contains all four colours.  Propagation from all 24 assignments on one seed
tetrahedron finds exactly 24 complete labelled colourings on each carrier:
the colouring is unique up to global renaming, with no later local branch.

All 120 left multiplications, 120 right multiplications and quaternion
conjugation preserve each declared colouring, with zero failures among 241
tests per carrier.

These are **DERIVED COMPUTATIONAL** combinatorial facts.

## 2. The slabs really are conforming

Using the declared representative order `0<1<2<3`, the standard staircase
product gives:

| carrier | pentachora | distinct four-faces | boundary four-faces | interior four-faces |
|---|---:|---:|---:|---:|
| `P(sd K)` | 57,600 | 158,400 | 28,800 | 129,600 |
| `P(Esd_2(sd K))` | 460,800 | 1,267,200 | 230,400 | 1,036,800 |

Every four-face has incidence one or two.  The incidence-one set is exactly
the bottom and top copies of the input closed spatial complex; there is no
side boundary.  All pentachora are distinct.

Thus the previous obstacle was not existence or conformity.  A compact,
globally compatible fine slab can be built.

## 3. Why existence does not select the slab

For a fixed proper four-colouring, each permutation of the colour classes is
a valid colour-consecutive order.  The complete census gives

```text
N_order = 24 distinct labelled slab complexes,
N_H4    = 24 spatially H4-invariant alternatives.
```

Time reversal sends an order

```text
(a,b,c,d) -> (d,c,b,a).
```

There are zero fixed orders and twelve orbits of size two.  Declaring a time
orientation therefore does not select one spatial colour order.

The audit encoded an order by the diagonal it places on each of the six
boundary squares `edge x I` of a tetrahedral prism.  These six choices form a
tournament on four colours.  Of all `2^6=64` tournaments, exactly 24 are
transitive, and they are exactly the tournaments of the 24 total orders.
Because the six diagonals are recoverable from the slab, the encoding is
injective: the 24 global schedules cannot be an artifact of a pentachoron
set comparison.

The audit also supplied two controls:

- a segment prism has exactly its two expected diagonals, exchanged by time
  reversal;
- a tournament containing a directed three-cycle is rejected as no total
  staircase order.

Finally, permutation parity divides the schedules as `12+12`.  Even adding a
spatial orientation could at best halve the family, not select one member.

This count is **DERIVED COMPUTATIONAL, ADVERSARIALLY CORROBORATED**.

## 4. The substantive judgement: what data are actually present?

It is essential to separate:

```text
D0  uncoloured carrier                              -> 24 schedules
D1  unordered partition into four colour classes   -> 24 schedules
D2  map into the ordered set 0<1<2<3               ->  1 schedule.
```

`D2` selects one slab tautologically.  The issue is whether its linear order
was already forced by the certified spatial geometry.  The pre-mission
carrier artifact certified rank subdivision but contained no temporal,
slab, colour-order or schedule datum.  The fine classes are modular
weighted-rank classes, not literal face dimensions.  Declaring their residue
representatives to carry the temporal order `0<1<2<3` therefore adds the
missing `D2` structure.

A defensible counterargument is that the standard edgewise construction is
performed in ordered rank coordinates, so its modular colour labels should
inherit the ordinary order of their representatives.  I do not know a
theorem that promotes that bookkeeping order to a geometrically forced
temporal diagonal choice.  The product literature treats a linear vertex
order as input and warns that it matters.  Consequently:

- under the frozen protocol, the carrier is **STRUCTURAL**;
- whether a stronger functorial notion should regard the rank-residue order
  as already part of the carrier is **OPEN**;
- it cannot be upgraded after seeing the census without a new protocol and a
  new selection axiom.

This is the hostile reading.  It prevents a convenient arithmetic convention
from being promoted silently into a physical tick.

## 5. Literature reconciliation

The post-result search reinforced rather than removed the distinction.

Joswig and Witte define the simplicial product using linear orders on the
vertices of both factors.  In Section 3 they state that the ordering is
crucial and exhibit different orderings giving pairwise non-isomorphic
product triangulations.  Proposition 3.2 shows that colour-consecutive orders
give a balanced product; it does not choose the order:

- M. Joswig and N. Witte, *Products of Foldable Triangulations*,
  arXiv:`math/0508180v3`.

The local edgewise carrier and its coloured-barycentric setting are standard:

- H. Edelsbrunner and D. R. Grayson, *Edgewise Subdivision of a Simplex*,
  DOI `10.1007/s004540010063`;
- C. A. Athanasiadis, *Edgewise subdivisions, local h-polynomials and
  excedances in the wreath product*, arXiv:`1310.0521`, DOI
  `10.1137/130939948`.

No external novelty is claimed for either construction.  The search did not
locate a primary source deriving a unique temporal staircase order from the
present modulo-four rank colouring.  Search absence is not proof, so external
novelty and the strongest functorial selection question remain **OPEN**.

## 6. Physical status and next boundary

- **DERIVED:** a compact globally conforming slab exists on both canonical
  spatial carriers.
- **DERIVED NEGATIVE under the preregistered selection rule:** the spatial
  carrier and `H4` symmetry do not select one of the 24 staircase schedules.
- **STRUCTURAL:** choosing `0<1<2<3`, choosing any other schedule, or averaging
  schedules without a new principle.
- **OPEN:** a schedule-free polytopal/cellular Regge action on tetrahedral
  prisms, or an independently derived temporal structure that supplies `D2`.
- **NOT RUN:** local lapse equations, anisotropic Hessian, dispersion and
  effective limiting speed on this slab family.

The route is not globally dead: the spatial carrier and local P1 dust weights
remain valid.  But the simplicial temporal-carrier branch stops here.  Running
one schedule would manufacture the tick whose derivation is the point of the
programme.

Only the two mission verifiers and static guards are run.  The full suite is
not run by explicit user instruction.

