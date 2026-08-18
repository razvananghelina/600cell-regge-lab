# A canonical four-dimensional slab exists, but it has no selected dynamics

Date: 2026-08-12

Preregistered protocol commit: `5c9fff8`

Registered verifier:
`reproducible/verify_gravity_time_slab_canonicity.py`

Machine-readable result:
`reproducible/gravity_time_slab_canonicity.json`

## Headline

There are three distinct results.

> **DERIVED CANONICAL CW SLAB.** The product `K x I` is a choice-free,
> `H4`-equivariant four-dimensional cellular carrier with two spatial
> boundary copies of the 600-cell boundary `K`.

> **DERIVED SIMPLICIAL CHOICE OBSTRUCTION.** No triangulation of that product
> using only its product vertices can preserve full `H4`. Even double
> transpositions in the rotational tetrahedron stabilizer exchange the two
> possible diagonals of an `edge x I` square. Barycentric subdivision repairs
> canonicity by adding all cell barycentres.

> **DERIVED CURRENT SLAB-DYNAMICS GAP.** The barycentric cylinder has exactly
> 115,200 maximal four-flags, equal to the state count of the existing robust
> walk, but none of all `8!` projection-preserving component identifications
> makes its literal stages into cylinder-chamber incidence moves.

The targeted verifier passes `23/23`. No full suite was run.

The cardinality equality and the likely factorized failure were noticed before
preregistration and disclosed in the protocol. The exhaustive incidence table,
not the number 115,200, is the result.

## 1. Cap, slice and slab are different objects

Let

```text
f(K)=(120,720,1200,600).
```

The cone has

```text
f(CK)=(121,840,1920,1800,600),     chi(CK)=1.
```

It is the canonical `H4`-symmetric four-ball obtained by joining the centre to
every boundary tetrahedron. Its only boundary is `K`. It can be read as a
Euclidean filling or state-preparation cap, but it is not a map between two
spatial slices.

The product CW complex instead has

```text
f(K x I)=(240,1560,3120,2400,600),     chi(K x I)=0.
```

Its 600 four-cells are tetrahedral prisms and its boundary contains two exact
copies of `K`. The 120 vertical edges form one `H4` orbit. Thus the geometry
provides a natural carrier for vertex lapse data.

This is the first honest place in the construction where “two consecutive
spatial frames” exist. It is still only topology. The interval has no selected
proper length, signature or action, and a vertical-edge variable has not been
shown to act as a Lagrange multiplier.

## 2. Why a minimal simplicial product is not canonical

The stabilizer of one tetrahedron has order 24 and induces the complete `S4`
on its vertices. The obstruction already occurs in the even subgroup: a
double transposition `(ij)(kl)` swaps the endpoints of the edge `(ij)`.

The product of that edge with `I` is a square. A simplicial triangulation using
only its four vertices must choose exactly one of its two diagonals, while the
double transposition exchanges them. Therefore no such choice is invariant.

Using an even element is load-bearing. Pairing odd spatial reflections with
time reversal cannot repair the obstruction, because these double
transpositions have positive orientation character and still exchange the
diagonals.

Hence a total vertex ordering or a local diagonal convention would break
`H4`. It may be useful numerically, but it would not be selected by the theory.

Barycentric subdivision is the canonical escape: every nonempty product cell
gets a vertex, so no diagonal is chosen. It has

```text
7,920 cell-barycentre vertices
115,200 maximal four-dimensional flag simplices
28,800 boundary flag tetrahedra across the two copies of sd(K).
```

The maximal count factors exactly as

```text
600 prisms * 24 spatial flags * 2 interval endpoints * 4 shuffles
=115,200.
```

The price is real: the boundary carrier is now `sd(K)`, not unrefined `K`, and
the metric data of the added vertices are not fixed by the face poset alone.

## 3. Exact product-chamber geometry

A maximal flag in one tetrahedral prism can be labelled by

```text
(spatial chamber, epsilon, r),
epsilon in {0,1},       r in {0,1,2,3}.
```

Here `epsilon` is the interval endpoint at which the flag starts, and `r` is
where the interval promotion is interleaved with the three spatial
promotions.

The eight labels form two length-four paths, joined only between their `r=0`
endpoints. A spatial chamber-colour move `s_i` lifts to one product-chamber
edge exactly when it preserves `(epsilon,r)` and `r != i`. At `r=i`, the
spatial face being changed occurs twice in the product flag, so two ranks
would change rather than one.

The verifier reconstructs the local face chains themselves and checks these
rules; they are not inferred from the robust walk.

## 4. Why 115,200 was a false bridge

The robust walk has

```text
14,400 spatial chambers * 8 components = 115,200 states.
```

The cylinder has

```text
14,400 spatial chambers * 8 product-flag labels = 115,200 chambers.
```

Cardinality therefore makes a factorized identification tempting. There are
exactly `8!=40,320` component-to-product-label bijections. Every one was
tested with identity/idling allowed.

The exact outcome is:

```text
stage S0 local for       960 / 40,320 labelings
stage S2 local for     2,880 / 40,320 labelings
stage S1 local for         0 / 40,320 labelings
macro one-edge local       0 / 40,320 labelings
all three stages           0 / 40,320 labelings.
```

The best `S1` and macro identifications satisfy six of eight transition rules.
The failure is structural: their cross-chamber transitions simultaneously
change the robust component, whereas a product-chamber edge crossing a
spatial facet preserves `(epsilon,r)`.

The macro one-edge test is a diagnostic, not the load-bearing negative: a
composition of three local stages need not itself be one edge. The decisive
fact is that the literal intermediate stage `S1` is nonlocal for every
projection-preserving identification.

This does not rule out a chamber-dependent relabelling. Such a relabelling
would alter the projection to the spatial chamber and needs a new geometric
selector; it cannot be introduced merely to save the equality of dimensions.

## 5. Pachner moves on the initial slice

The exact initial legal-move census is

| move | legal sites | `H4` orbits |
|---|---:|---:|
| `1-4` | 600 | 1 |
| `2-3` | 1,200 | 1 |
| `3-2` | 0 | 0 |
| `4-1` | 0 | 0 |

All edges initially have five incident tetrahedra, so no `3-2` move is legal;
all vertices have twenty incident tetrahedra, so no `4-1` move is legal.

The complete `1-4` orbit is compatible and gives a canonical stellar
refinement with

```text
f=(720,3120,4800,2400).
```

It is a refined spatial slice rather than a translation between equal
carriers.

Every legal `2-3` move conflicts with six others because their two-tetrahedron
supports overlap. Since all 1,200 sites form one transitive orbit, the only
full-`H4` invariant subsets are empty and full, and the full set is not a
parallel layer. Thus no nonempty `H4`-invariant one-layer `2-3` evolution
exists.

This does not exclude a longer orbit-complete schedule or a quantum sum over
moves. Either would need a new ordering/amplitude principle.

The physical interpretation is external and conditional. Canonical
simplicial gravity treats simplex gluing actions as generators of discrete
evolution [Dittrich--Hoehn](https://arxiv.org/abs/1108.1974). On a flat
linearized background, `1-4` creates lapse/shift gauge variables while `2-3`
creates a lattice graviton [Hoehn](https://arxiv.org/abs/1411.5672). The
present slice is not such a flat spacetime background. In curved finite Regge
calculus exact gauge symmetries can be broken into pseudo-constraints
[Bahr--Dittrich](https://arxiv.org/abs/0905.1670), so none of those physical
roles may simply be copied here.

## 6. Status ledger

| Claim | Status |
|---|---|
| `CK` is a canonical four-dimensional filling | **DERIVED** |
| `CK` is a two-slice time evolution | **REFUTED** |
| `K x I` is a canonical two-boundary CW carrier | **DERIVED** |
| Its 120 vertical edges form a canonical lapse carrier | **DERIVED CARRIER ONLY** |
| A vertical edge is already a lapse multiplier | **OPEN / NOT DERIVED** |
| A no-new-vertex simplicial product preserves `H4` | **REFUTED** |
| `sd(K x I)` is a canonical simplicial carrier | **DERIVED** |
| Its four-chamber count is 115,200 | **DERIVED** |
| The robust 115,200-state walk is its factorized chamber walk | **REFUTED** |
| The complete `1-4` orbit is a canonical stellar refinement | **DERIVED** |
| A nonempty invariant parallel `2-3` layer exists | **REFUTED** |
| A longer selected Pachner schedule exists | **OPEN** |
| A Lorentzian metric/action/constraint algebra exists on the slab | **OPEN** |
| The 47-dimensional kinetic family has been reduced | **OPEN** |
| A fourth physical dimension, `c`, `G` or Planck time follows | **NOT CLAIMED** |

## 7. Physical consequence and next gate

The phrase “we only had a spatial frame” can now be sharpened:

> We have constructed a canonical carrier containing two frames and one
> combinatorial interval, but not the law that transports the metric from one
> frame to the other.

The next admissible gate is metric and variational, not another state-count
match. On `K x I` or its canonical barycentric subdivision it must specify,
before examining the 150 negative spatial modes:

1. the independent spatial, vertical and cross-edge metric data;
2. Euclidean versus Lorentzian simplex admissibility;
3. one local action on every four-cell;
4. its discrete Legendre transform;
5. lapse/shift pre- and post-constraints and their class;
6. whether the resulting kinetic block selects a point in the 47-dimensional
   `H4` family.

The cleanest candidate supplied by established discrete gravity is a Regge
action used as Hamilton's principal function. It is an **external candidate**,
not yet a derived theory choice. It must be preregistered together with every
metric/triangulation choice before any conformal-mode comparison.

The subsequent local audit `gravity_tent_move_regge_result.md` implements
that candidate on the canonical icosahedral vertex tent. Under Euclidean,
zero-volume-coefficient and static-equilateral hypotheses it derives
`t/a=phi^-1` exactly from fivefold hinge incidence. Releasing the final
spatial length reveals a continuous flat family, so the result is a
conditional static return ratio rather than a selected physical tick.

## 8. Reproduction history

The first targeted run passed `22/22`. Hostile review then noticed that an odd
spatial reflection might be paired with reversal of `I`. The triangulation
obstruction was strengthened to use the three even double transpositions in
the rotational tetrahedron stabilizer; a new exact check identifies that
stabilizer as `A4`, so time reversal cannot repair the diagonal choice. The
exact counts and all `8!` comparison results were unchanged. The final
targeted run passes `23/23`.

No full suite and no PDF build were run.
