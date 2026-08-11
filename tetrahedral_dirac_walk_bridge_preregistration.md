# Preregistration: ordered barycentric chambers versus the tetrahedral Dirac walk

Date: 2026-08-11

## External construction being tested

Nzongani, Eon, Márquez-Martín, Pérez, Di Molfetta and Arrighi,
*Dirac quantum walk on tetrahedra* (2024), construct a local causal unitary
walk on a flat tessellation by left- and right-handed tetrahedral
orthoschemes.  Their carrier has four facet amplitudes per tetrahedron.  The
scheme requires:

1. four ordered neighbours, one across each labelled facet;
2. a left/right handedness label;
3. a prescribed local permutation shift;
4. spin coins chosen so that the continuum expansion is the Dirac equation.

The paper proves the (3+1)-dimensional Dirac limit for its flat periodic
tessellation.  It explicitly leaves propagation on a curved triangulated
manifold as an open question.  Therefore its continuum theorem will **not**
be transferred to the 600-cell.

## Target-blind geometric question

Take the boundary complex of the regular 600-cell and its first full
barycentric subdivision.  Its top simplices are complete flags

\[
v\subset e\subset f\subset t.
\]

Before constructing a spectrum or comparing with a continuum equation, test
whether these flags canonically supply the discrete input required by the
published walk:

- a four-regular dual chamber graph;
- an intrinsic edge colour (0,1,2,3), recording which flag rank changes;
- fixed-point-free involutions (s_0,s_1,s_2,s_3);
- the Coxeter relations of type (H_4), with adjacent orders (3,3,5) and
  nonadjacent order (2);
- a balanced bipartition interpreted only as chamber orientation;
- bijectivity of the published two-stage shift after identifying its four
  facet labels with the four intrinsic rank colours.

No Standard-Model number, particle mass, speed, Planck scale or continuum
dispersion is used.

## Construction frozen before execution

1. Rebuild the `120/720/1200/600` 600-cell complex from the shared canonical
   vertex construction.
2. Enumerate all 24 ordered flags inside every coarse tetrahedron, giving
   14,400 chambers.
3. Define (s_0,s_1,s_2) by swapping adjacent entries of the flag ordering.
4. Define (s_3) by retaining the ordered boundary triangle and replacing
   the parent tetrahedron by the unique coarse tetrahedron across it.
5. Build the undirected, rank-coloured chamber graph from these four maps.
6. Compute a bipartition from the graph itself; do not prescribe handedness
   from a desired walk.
7. On the carrier `(chamber, facet component)` of dimension 57,600, implement
   exactly the two permutation stages (S_B,S_G) in Eqs. (2)--(3) of the
   paper, using colour 2 across one handedness and colour 3 across the other.
   Check the maps as integer permutations rather than floating matrices.

## Decision boundaries

- **DERIVED STRUCTURAL BRIDGE:** all coloured chamber axioms and the
  published shift bijectivity hold.  Then the 600-cell refinement supplies a
  canonical ordered/chiral carrier on which that shift can be run.
- **DERIVED NEGATIVE:** a colour map is not an involution, the Coxeter
  relations fail, the graph is not balanced bipartite, or the shift is not a
  permutation.  Then this direct bridge is closed.

Passing is deliberately not called a Dirac derivation.  The local spin coin
was designed in the paper to recover the Dirac equation, and convergence was
proved only for a flat periodic orthoscheme lattice.  On the closed curved
600-cell carrier, both selection of the complete dynamics and a refinement
limit remain **OPEN**.

## Hostile canonicity audit

Rank colours are intrinsic because the four flag elements have different
dimensions.  The global exchange of the two bipartition labels is an
orientation convention and must be reported.

The **shift rule itself is not derived from 600-cell geometry**.  It is an
external, published candidate transplanted onto a carrier that the geometry
selects.  Even a passing result is therefore STRUCTURAL, not a new physical
law.  A later acceptance gate must show either that the walk is selected by
the theory's own operator/axioms or that all remaining choices are gauge
equivalent.
