# Prior-art gate: balanced temporal slab on the projected rank-edgewise carrier

Date: 2026-08-19

## Exact question and hypotheses

Let

```text
K0 = P(sd K_600),
K1 = P(Esd_2(sd K_600)),       P(x)=x/||x||,
```

be the two already certified closed spatial carriers.  A vertex of `K0` is a
nonempty face of the 600-cell and has its intrinsic rank

```text
r(v)=dim(face) in {0,1,2,3}.
```

A vertex of `K1` is either an old vertex `(v,v)` or the midpoint `(v,w)` of
an edge of `K0`.  The candidate colouring is

```text
c(v,w) = r(v)+r(w) mod 4.
```

The candidate temporal carrier is the staircase triangulation of `Ki x I`
obtained by ordering the four colours as `0<1<2<3` in every spatial
tetrahedron.  Time orientation is stated explicitly: the two interval
endpoints are ordered `past<future`.

The exact gate is:

1. is `c` a proper global four-colouring, with every tetrahedron containing
   all four colours;
2. is the colouring forced by the abstract carrier up to one global `S4`
   relabelling, or are there further colourings;
3. does the rank formula fix that global relabelling without vertex IDs,
   metric comparisons or a fitted schedule;
4. does the induced staircase rule give one conforming triangulation of
   `Ki x I`, with exactly the two copies of `Ki` as boundary;
5. is it preserved by the declared spatial `H4` action;
6. under reversal of time orientation, is the slab identical or a distinct
   legitimate staircase triangulation.

No Regge action, lapse equation, Hessian, continuum coefficient or physical
target belongs to this gate.  The carrier must be settled first.

## Known construction

The construction is not a new triangulation principle.

Edelsbrunner and Grayson give the standard edgewise subdivision of a simplex,
including its colour-scheme description:

- H. Edelsbrunner and D. R. Grayson, *Edgewise Subdivision of a Simplex*,
  DOI `10.1007/s004540010063`.

Athanasiadis studies the edgewise subdivision of a barycentric subdivision,
also called a coloured barycentric subdivision:

- C. A. Athanasiadis, *Edgewise subdivisions, local h-polynomials and
  excedances in the wreath product*, arXiv:`1310.0521`, DOI
  `10.1137/130939948`.

The directly relevant product theorem is:

- M. Joswig and N. Witte, *Products of Foldable Triangulations*,
  arXiv:`math/0508180v3`.

Their Section 3 defines the simplicial/staircase product of two complexes.
It proves that it triangulates the product, but also states explicitly that
the linear orderings of the vertices are crucial and may produce pairwise
non-isomorphic triangulations.  Proposition 3.2 proves that colour-consecutive
orders on foldable (balanced) factors give a foldable product.  Equation (3)
shows that, on each product cell, the lifting perturbation can depend only on
the colour classes.  The same paper notes that a foldable pure complex with
connected dual graph has its minimal colouring uniquely up to renaming the
colours.

Thus the literature supplies the existence and conformity theorem once a
global ordered colouring is part of the data.  It does **not** select that
ordered colouring for the present 600-cell carrier.

## Repository controls and excluded branches

The repository already has:

- a canonical rank order on every barycentric chamber and a certified first
  rank-edgewise refinement;
- a barycentric product slab which is not a common refinement of this carrier;
- a universal staircase overlay which is compatible but very large and only
  **STRUCTURAL**;
- a different five-colour schedule branch whose remaining chiral `Z2`
  ambiguity is not selected by `gamma` or `J`.

The present four-colouring comes from face ranks on the certified spatial
carrier.  It is not the five-colour coset cover, so the old no-go is not
reused as a theorem here.  It is retained as a warning that a named schedule
can hide a discrete choice.

## Preliminary observation and provenance limitation

Before this formal gate, a scratch calculation on one rank-edgewise chamber
found that `c(v,w)=r(v)+r(w) mod 4` colours all eight children properly and
that the local ten-vertex graph has 24 labelled proper colourings.  This is
disclosed rather than presented as preregistered evidence.  The registered
test must reconstruct the complete carriers, enumerate their colourings and
test the global product topology.

## Classification before execution

- **KNOWN:** staircase products triangulate `|K| x I`; the result depends on
  vertex order; colour-consecutive orders give a balanced product.
- **DERIVED INPUT:** face rank and the rank-edgewise carrier are invariant
  under the certified 600-cell symmetry action.
- **CONTROL:** `K0` must recover its ordinary face-rank colouring and both
  slab boundaries must reproduce the input spatial complex.
- **OPEN:** whether the fine formula is proper globally and unique up to
  `S4`; whether it selects a single oriented slab.
- **OPEN:** external novelty of using this exact carrier in an irregular
  Regge--dust evolution.  Search absence cannot establish novelty.

## Decision rule

The route advances only if the rank data select a globally conforming slab
once a time orientation is declared.  If an additional vertex order,
colour permutation, diagonal schedule or metric comparison changes the slab,
the construction is labelled **STRUCTURAL** and cannot yet support a derived
local lapse.  Time reversal may legitimately produce the oppositely oriented
slab, but that dependence must be counted and stated rather than hidden.

