# Prior-art gate: chromatic degree as a staircase-parity selector

Date: 2026-08-17

Status: written before computing the degree of the fixed 600-cell colouring or
comparing its sign across the 120 staircase orders.

## Exact object

Fix the already certified oriented 600-cell boundary complex `K` and its five
24-vertex colour classes.  Since every tetrahedron contains four distinct
colours, the colouring defines a simplicial map

```text
f : K -> boundary(Delta^4),
```

where a tetrahedron missing colour `i` maps to the facet of the standard
four-simplex missing vertex `i`.

Orient `K` with the frozen determinant chain.  Orient `Delta^4` by the vertex
order `(0,1,2,3,4)` and its boundary by

```text
boundary[0,1,2,3,4]
  = sum_i (-1)^i [0,...,omit i,...,4].
```

The integer degree of `f` is obtained by pushing the complete signed
tetrahedral fundamental chain forward.  It is not a numerical fit and does not
use the Regge action.

For a staircase order `sigma`, relabel each colour by its rank in `sigma` and
compute the degree again.  The mathematical expectation is alternation by
`sign(sigma)`, but the base degree may be zero; neither outcome may be assumed.

## Primary prior art

- Fisk proves that the 600-cell has exactly ten five-colourings and describes
  them through left and right binary-tetrahedral cosets: [Coloring the 600
  Cell](https://arxiv.org/abs/0802.2533).
- The degree of a colouring map to the boundary of a simplex is a standard
  combinatorial degree and is alternating in the colour order: [Toric residue
  and combinatorial degree](https://arxiv.org/abs/math/0309409), especially
  Definition 1.11.
- Staircase triangulations of a simplex prism correspond to total orders and
  adjacent transpositions: [Non-connected toric Hilbert
  schemes](https://arxiv.org/abs/math/0204044), Propositions 1.2--1.3.

The literature establishes the mathematical invariant.  No located source
states that a positive chromatic degree is a physical causality condition or a
Regge time-order selector.

## KNOWN / CONTROL / OPEN

### KNOWN

- All 120 staircase orders are distinct, coherently oriented triangulations
  with the same relative product chain.
- The exact setwise `H4` action on the fixed cover induces `A5`, producing one
  even and one odd schedule orbit.
- Relabelling the five target vertices by an odd permutation reverses the
  target orientation and hence the sign of any nonzero degree.

### CONTROL

- Compute the pushforward coefficient separately on each of the five target
  facets; after multiplying by the target boundary sign, all five integers
  must agree.
- Independently compute the degree as a signed preimage count of one target
  facet.
- Repeat over all 120 relabellings and require exactly
  `degree_sigma = sign(sigma)*degree_identity`.
- Require invariance under all 60 induced `A5` cover permutations.
- Load no Regge, nonlinear, continuum or desired-parity target.

### OPEN

- Whether the fixed colouring degree is zero or nonzero.
- If nonzero, its integer magnitude.
- Whether the two schedule orbits are exactly its two signs.
- Whether any already-derived physical axiom requires one sign rather than the
  other.

## Decision boundary

- `CHROMATIC_DEGREE_ZERO`: the map supplies no orientation line.
- `CHROMATIC_ORIENTATION_LINE_DERIVED`: the degree is resolved nonzero, all
  controls pass and its sign separates the two schedule orbits.
- `OPEN_CONTROL_FAILURE`: the facet degrees disagree or any exact control
  fails.

Even `CHROMATIC_ORIENTATION_LINE_DERIVED` is only a **STRUCTURAL** selector.
Call it a physical schedule selection only if a separately stated, already
derived dynamical or causal axiom requires positive degree.  Choosing the sign
after seeing which Regge schedule is favourable is forbidden fitting.

External novelty remains **OPEN**.
