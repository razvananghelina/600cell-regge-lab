# Prior-art gate: anisotropic rigidity of a tetrahedral spacetime frustum

Date: 2026-08-19

## Exact question and complete hypotheses

Consider one flat Lorentzian four-dimensional block with the combinatorics

```text
tetrahedron x interval.
```

Its boundary graph has four bottom vertices, four top vertices, six bottom
edges, six top edges and four corresponding struts.  The proposed
schedule-free cellular extension treats those 16 squared lengths as the
complete local metric data and asks for the Regge hinge areas and dihedral
angles needed by an anisotropic Hessian.

The question is:

> Do the 16 cellular graph lengths determine the flat four-dimensional
> tetrahedral frustum locally, modulo Lorentz isometries, or are additional
> internal shape data required?

This question precedes dust, curvature, the 600-cell group action and every
physical spectrum.  A local failure kills the proposed anisotropic cellular
Hessian even if the homogeneous frustum action is valid.

## Disclosed dimension count

Eight generic vertices in four dimensions have

```text
8*4 - dim ISO(3,1) = 32 - 10 = 22
```

local geometric degrees of freedom modulo translations and Lorentz
transformations.  The cellular graph supplies only

```text
6 bottom + 6 top + 4 struts = 16
```

squared lengths.  The generic deficit is therefore at least six.  This
count is **DERIVED**, but its applicability at the symmetric frustum used by
the repository is still to be tested exactly; symmetry can increase rather
than decrease infinitesimal degeneracy.

Each staircase triangulation chooses one diagonal on each of the six
quadrilateral faces `edge x interval`.  It therefore supplies exactly six
additional cross-edge lengths.  The known 24 colour orders give 24 such
six-diagonal completions.  The target-disclosed prediction is:

```text
cellular graph rigidity rank          16,
non-isometric infinitesimal flexes      6,
every staircase completion rank        22.
```

No rigidity matrix has been evaluated while writing this gate.

## What primary literature establishes

The Collins--Williams/polytopal programme uses frusta in symmetry-reduced
cosmological models.  Tsuda and Fujiwara vary common edge lengths and struts
to obtain a Hamiltonian constraint and an evolution equation; their
construction does not claim that arbitrary anisotropic tetrahedral-frustum
geometry is encoded by the 16 graph lengths:

- R. Tsuda and T. Fujiwara, *[Higher Dimensional Polytopal Universe in Regge
  Calculus](https://arxiv.org/abs/2109.01075)*;
- R. Tsuda and T. Fujiwara, *[Oscillating 4-Polytopal Universe in Regge
  Calculus](https://arxiv.org/abs/2011.04120)*.

Their earlier polytopal treatment notes that a non-simplicial flat block can
be fully triangulated by adding hinges with vanishing deficit.  This shows
how to evaluate a specified flat block; it does not show that the unaugmented
boundary graph selects its missing shape:

- R. Tsuda and T. Fujiwara, *[Expanding polyhedral universe in Regge
  calculus](https://academic.oup.com/ptep/article/2017/7/073E01/4056197)*.

Standard length Regge calculus instead fixes a triangulation and uses its
edge lengths as the metric variables.  Dittrich and Steinhaus show that
triangulation independence in four-dimensional linearized Regge calculus is
not automatic; changing a triangulation is not generally an innocuous
classical relabelling:

- B. Dittrich and S. Steinhaus, *[Path integral measure and triangulation
  independence in discrete gravity](https://arxiv.org/abs/1110.6866)*.

Barrett's first-order Regge formulation makes the complementary point that
extra angle/area variables can be introduced, but then they are additional
geometric variables with compatibility conditions rather than consequences
of an insufficient edge set:

- J. W. Barrett, *[First order Regge
  calculus](https://arxiv.org/abs/hep-th/9404124)*.

No primary source located in the 2026-08-19 search asserts local rigidity of
the 16-edge tetrahedral-prism graph in four dimensions.  Search absence is
not a proof; the exact local rank calculation below is the evidence.

## Repository context

- **DERIVED:** the homogeneous cellular action is well defined on the
  projected refinement tower because each block is restricted to common
  lower/upper scale and common lapse.
- **DERIVED NEGATIVE:** the canonical projected spatial carrier admits 24
  equally symmetric staircase time orders; spatial `H4` plus time orientation
  selects none.
- **OPEN AT THIS GATE:** whether the cellular block avoids that ambiguity by
  being locally rigid without the six selected diagonals.

The homothetic ansatz must not be used as evidence for anisotropic rigidity:
it fixes relative orientation and shape before variation, precisely the data
under investigation.

## Acceptance and kill boundary

- If the 16-edge graph has no non-isometric flex and the count above fails,
  the schedule-free cellular route survives and its metric data must be
  characterized.
- If it has six or more flexes and each staircase completion removes them,
  the proposed cellular anisotropic action is underdetermined.  Running a
  Hessian then requires either selecting a triangulation or adding new
  area/angle/shape variables.
- Selecting one of 24 completions remains **STRUCTURAL**.
- Uniformly averaging their actions is a different ensemble theory and is
  not authorized by this classical rigidity test.

Even a clean negative does not refute Regge calculus.  It closes only the
claim that the repository's current schedule-free cellular boundary data
already determine a refined anisotropic classical tick.

