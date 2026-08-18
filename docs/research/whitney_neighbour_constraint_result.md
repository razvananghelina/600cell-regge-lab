# Uniformly local Whitney spectral constraints survive refinement

Date: 2026-08-11

Preregistration commit: `a819a52`

Targeted verifier:
`reproducible/verify_whitney_neighbour_constraints.py`

Targeted result: **12/12 PASS**.  The verifier is registered.  The full suite
was not run by explicit user request.

## Result

Replace complete-star averaging in the Whitney KKT pencil by pairwise
constraints between tetrahedron copies that share a triangle containing the
copied simplex.

For a (p)-simplex (s), the occurrence nodes are ​((t,s)) with
(s\subset t).  The canonical constraint graph joins

\[
(t,s)\longleftrightarrow(t',s)

\]

exactly when (t\cap t') is a triangle containing (s).  Every constraint
row is an unweighted signed difference `(1,-1)`; there are no anchors,
spanning trees, weights or fitted coefficients.

The complete exact result is:

| level | local cochains | constraint rows | exact rank | redundant multiplier gauge | maximum copy degrees |
|---|---:|---:|---:|---:|---|
| base 600-cell | 9,000 | 8,400 | 6,360 | 2,040 | `(3,2,1,0)` |
| first barycentric | 216,000 | 201,600 | 153,120 | 48,480 | `(3,2,1,0)` |

At both levels the verifier proves:

\[
CJ=0,
\qquad
\ker C=\operatorname{im}J,
\qquad
\operatorname{rank}C=N_{\rm loc}-N_{\rm global}.
\]

> **DERIVED UNIFORMLY LOCAL SPECTRAL CONSTRAINT:** exact conforming assembly
> can be imposed using degree-at-most-three neighbour equations, independently
> of the growing simplex-star size.

Together with commit `dc7f3df`, this gives an exact bounded-degree KKT pencil
for the Whitney generalized spectrum without applying the assembled mass
inverse.

## Exact kernel proof

Every closed triangulated 3-manifold triangle has exactly two parent
tetrahedra.  For each shared triangle, the construction adds constraints for
its three vertices, three edges and the triangle itself, hence seven rows.
The exact totals are therefore

\[
7\times1200=8400

\]

and

\[
7\times28800=201600.
\]

Constraint edges never mix two global simplices.  Inside the copies of one
simplex, they form the dual adjacency graph of its tetrahedral star.  The
verifier enumerates every such graph and proves it connected at both levels.
A signed graph incidence matrix has rank `nodes - components`, so connectivity
gives the kernel and ranks exactly without a numerical rank tolerance.

The rank totals agree with

\[
9000-2640=6360

\]

and

\[
216000-62880=153120.
\]

## Why the degree remains bounded

A tetrahedron copy containing a fixed (p)-face has exactly (3-p)
codimension-one faces which also contain that (p)-face.  Across each such
face lies one neighbouring tetrahedron.  Thus the occurrence-node degrees
are at most

\[
(3,2,1,0)

\]

for form degrees (0,1,2,3).  Both enumerated levels saturate these values.

This extends beyond the two computations under complete hypotheses:

> **DERIVED GENERAL STATEMENT:** on a closed combinatorial 3-manifold whose
> vertex and edge links are connected, with exactly two tetrahedra at every
> triangle, the face-neighbour constraints have kernel equal to conforming
> assembly and degree at most three.  Barycentric subdivision preserves these
> link properties, so the bound persists at every refinement level.

The computations certify that the actual 600-cell and its first subdivision
satisfy the hypotheses; the general step is the graph-incidence argument
above.

## Spectral consequence

Because

\[
\ker C=\operatorname{im}J,

\]

finite-dimensional orthogonality gives

\[
\operatorname{ran}C^T=(\operatorname{im}J)^\perp.

\]

Therefore

\[
\begin{pmatrix}
A_{\rm loc}-zM_{\rm loc}&C^T\\
C&0
\end{pmatrix}

\]

has exactly the assembled Whitney physical finite spectrum.  The weak blocks
are tetrahedron-local and the constraints cross only one shared triangle.
Unlike the inverse polynomial, neither local degree nor algebraic stencil
depth grows under refinement.

This resolves the apparent contradiction:

- the **strong reduced generator** (M^{-1}A) becomes increasingly complex;
- the equivalent **weak constrained spectral pencil** remains uniformly
  local.

The inverse complexity was a consequence of eliminating local constraints,
not evidence that the underlying spectral relation itself must be nonlocal.

## Gauge redundancy and dynamical boundary

Keeping every incidence-selected neighbour row is canonical but redundant.
Cycles in the occurrence graphs give multiplier gauge dimensions 2,040 and
48,480.  Selecting a spanning tree would remove redundancy but introduce an
arbitrary noncanonical choice.

The multiplier block still has zero kinetic metric.  Consequently this is a
descriptor/constrained system, not an ordinary positive-metric Schrödinger
Hamiltonian.  The result therefore does not yet supply the physical tick.

It does suggest a sharper direction: treat copy equality as a genuine local
constraint system and perform a Dirac--Bergmann/gauge analysis rather than
eliminating the constraints into (M^{-1}).  Constrained local field theories
can have local equations even when their reduced instantaneous generators
look nonlocal.  Whether that happens here is **OPEN**, not assumed.

## Status ledger

- **DERIVED:** canonical neighbour constraint row counts at two levels.
- **DERIVED:** exact kernels, ranks and gauge redundancies.
- **DERIVED:** uniform maximum degree `(3,2,1,0)` under first refinement.
- **DERIVED GENERAL (stated hypotheses):** degree-at-most-three locality on
  every closed combinatorial 3-manifold with connected vertex/edge links and
  two tetrahedra per triangle.
- **DERIVED:** exact uniformly local representation of the Whitney spectrum.
- **DERIVED NEGATIVE:** the present singular descriptor is not an ordinary
  unitary tick.
- **STRUCTURAL:** interpreting constraint redundancy as a gauge-like sector.
- **OPEN:** Dirac--Bergmann classification and positive physical Hilbert
  evolution.
- **NOT CLAIMED:** Lorentzian time, inertia, mass or (c).

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_whitney_neighbour_constraints.py
```

Expected result: `12/12`.
