# Canonical barycentric product slab is not a staircase common refinement

Date: 2026-08-17

Only the targeted combinatorial verifier was run.  No gravity action and no
full suite were run.

## 1. Provenance

- prior-art gate: `dec110d`
- preregistered protocol: `f53857c`
- verifier registered before enumeration: `6cbbf5f`
- audited spatial source SHA-256:
  `a819ae9d472317d456cf7d67f588b31586b7aed2400a2540b8352f4661b39d45`
- result artifact:
  `reproducible/gravity_600cell_barycentric_product_slab.json`
- result SHA-256:
  `c35cf3ccca1fdcb4e144d3a0bea10c6955865b5660ff33abab3232e0ebf2315b`

The verifier passed `12/12` controls using exact incidence and rational affine
containment.

## 2. The canonical carrier that does exist

Let `K` be the 600-cell boundary and `X=K x I` its regular CW cylinder.  The
order complex of the product face poset is a canonical simplicial
triangulation of `S^3 x I`.  The exact f-vector is

```text
f(sd(K x I)) = (7920, 90480, 269760, 302400, 115200),
Euler characteristic = 0.
```

Its two induced boundary components each have

```text
(2640,17040,28800,14400),
Euler characteristic = 0.
```

The 2640 boundary vertices are precisely the nonempty cells of `K`, matching
the cardinality and incidence interpretation of the repository's oriented
cochain/Kähler--Dirac carrier.  Interval reflection is an exact simplicial
involution swapping the two boundaries.  Local `S4` vertex permutations and
time reflection preserve all product cells and all maximal chains.

These are **DERIVED COMBINATORIAL** facts.  They establish a schedule-free
four-dimensional carrier, not an action or dynamics.

## 3. The common-refinement test

A tetrahedral prism has 45 nonempty product cells and its barycentric order
complex has exactly 192 four-simplices.  Each was tested against the four
standard staircase simplices for every one of the 24 vertex orders.

The result is identical for all 24 orders:

```text
contained in one staircase simplex:  32 / 192
mixed across staircase simplices:    160 / 192
```

No barycentric four-simplex had multiple containing staircase simplices.  For
the canonical order, the 32 contained simplices split `16+16` between the two
middle staircase simplices; none lies wholly in either extreme simplex.  The
product-cell barycentre itself lies on the interface of the middle two but
outside both extreme staircase simplices, explaining why most barycentric
flags cross staircase diagonals.

Therefore the preregistered outcome is

```text
BARYCENTRIC_PRODUCT_NOT_A_STAIRCASE_COMMON_REFINEMENT.
```

## 4. Interpretation

- **DERIVED NEGATIVE:** the naive barycentric subdivision of the unsplit
  product CW complex cannot be used as a direct common refinement of even and
  odd staircase actions.
- **DERIVED:** the failure is universal over all 24 local vertex orders, not a
  peculiarity of the selected five-colouring.
- **STRUCTURAL POSITIVE:** it is nevertheless a canonical, time-reflection
  symmetric triangulation of the cylinder with the theory's 2640-cell carrier
  on each boundary.
- **NOT DERIVED:** a Lorentzian metric, nondegenerate four-simplex realization,
  dust action, canonical map or continuum dynamics on this carrier.

Replacing both schedules by this third triangulation would be a new
discretization choice, not coarse-graining evidence.  Calling it a perfect
action would be false.

## 5. Post-result prior art and the stronger construction

Adiprasito and Pak prove that any two PL-homeomorphic triangulations have a
common stellar subdivision: [*All triangulations have a common stellar
subdivision*](https://arxiv.org/abs/2404.05930).  This guarantees existence in
principle but does not provide a physically selected minimal refinement here.

Santos reviews geometric bistellar flips and secondary polytopes in
[*Geometric bistellar flips. The setting, the context and a
construction*](https://arxiv.org/abs/math/0601746), and studies triangulations
of products of simplices via the Cayley trick in
[*The Cayley trick and triangulations of products of
simplices*](https://arxiv.org/abs/math/0312069).

For the local prism, all staircase internal facets lie on the rational
hyperplanes

```text
t = sum_(i in A) lambda_i,
```

for nontrivial subsets `A` of the four tetrahedron vertices.  Their complete
arrangement defines the polyhedral overlay of all 24 staircases.  Barycentric
subdivision of that overlay, rather than of the unsplit prism, is a canonical
`S4 x C2`-invariant universal common refinement.

That statement identifies the next exact combinatorial object; its cell and
simplex census has not yet been computed here.

## 6. Status ledger

| Claim | Status |
|---|---|
| Product order complex triangulates `S^3 x I` | **KNOWN / independently reconstructed** |
| Printed global and boundary f-vectors | **DERIVED COMBINATORIAL** |
| Full local `S4 x C2` functoriality | **DERIVED COMBINATORIAL** |
| Naive product barycentric complex refines a staircase | **DERIVED NEGATIVE, 160/192 mixed** |
| It is a common refinement of even/odd schedules | **REFUTED** |
| Universal staircase-hyperplane overlay is finite and canonical | **STRUCTURAL** |
| Overlay census and global gluing | **OPEN** |
| Effective/perfect action on any refined carrier | **OPEN** |
| Unique nonlinear physical tick | **OPEN** |

## 7. Consequence

The cheap common-refinement proposal is closed, but it exposes the exact
missing geometry rather than returning to arbitrary scheduling: construct the
overlay of all staircase hyperplanes, then subdivide that overlay
functorially.  Only after that carrier is certified would it be honest to ask
for a coarse-grained Lorentzian Regge--dust action.
