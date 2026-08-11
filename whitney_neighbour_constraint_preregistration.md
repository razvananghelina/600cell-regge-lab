# Preregistration: uniformly local neighbour constraints for the Whitney pencil

Date: 2026-08-11

## Starting point and defect to repair

Commit `dc7f3df` proves that the exact Whitney generalized spectrum is the
finite spectrum of a local element KKT pencil with the equality constraint

\[
Q=I-J(J^TJ)^{-1}J^T.
\]

Although canonical, (Q) averages every tetrahedron copy in one simplex
star.  Its largest base block has size 20, while an original vertex belongs
to 120 child tetrahedra after one barycentric subdivision.  Calling this
uniform microscopic locality would therefore be false.

Replace (Q) by a constraint incidence matrix selected only by adjacent
top-dimensional cells.

## Frozen canonical graph

For every global (p)-simplex (s), its occurrence nodes are pairs

\[
(t,s),\qquad s\subset t,
\]

where (t) is a tetrahedron.  Join two occurrence nodes if and only if their
tetrahedra share a triangle (f) containing (s):

\[
t\cap t'=f,qquad s\subseteq f.
\]

Orient the constraint edge by the deterministic increasing tetrahedron index
and assign row ​((+1,-1)).  This orientation changes no kernel and is only a
storage convention; the undirected edge set is fixed by incidence.

Let (C_p) be the resulting signed graph incidence matrix.  It contains no
anchor, spanning-tree choice, weight or fitted coefficient.

On a closed tetrahedral 3-manifold, one top-cell copy has exactly

\[
3-p
\]

codimension-one faces containing its (p)-face.  Therefore the frozen local
degree bounds are `(3,2,1,0)` for (p=0,1,2,3).

## Two exact carriers

Construct and test both:

1. the base 600-cell boundary with f-vector `(120,720,1200,600)` and 9,000
   element-local cochains;
2. its complete first barycentric subdivision with f-vector
   `(2640,17040,28800,14400)` and 216,000 element-local cochains.

Every triangle must belong to exactly two tetrahedra at both levels.  Generate
constraints once per shared triangle and once per contained vertex, edge and
triangle.

The expected raw row counts, derived before enumeration, are

\[
7f_2=8400
\]

at the base and

\[
7f'_2=201600
\]

after refinement.

## Exact kernel and rank gates

For every level and degree:

1. every constraint row has exactly one `+1` and one `-1`;
2. both columns represent copies of the same global simplex;
3. (C_pJ_p=0) exactly;
4. the occurrence graph of every global simplex is connected;
5. distinct global simplices are never joined;
6. the maximum occurrence-node degree is exactly `(3,2,1,0)` whenever that
   degree has duplicated copies.

For a graph incidence matrix,

\[
\operatorname{rank}C_p
=N_p^{\rm loc}-\#\text{components}.
\]

Connectivity and nonmixing must therefore certify

\[
\ker C_p=\operatorname{im}J_p
\]

and

\[
\operatorname{rank}C_p
=N_p^{\rm loc}-N_p^{\rm global}
\]

without numerical rank tolerances.

Record redundant-row gauge dimensions

\[
\#\operatorname{rows}C_p-\operatorname{rank}C_p.
\]

Redundancy is allowed because retaining all incidence-selected neighbour rows
preserves canonicity; choosing a nonredundant spanning tree would reintroduce
an arbitrary choice.

## KKT consequence and boundary

If the gates pass, then

\[
\ker C=\operatorname{im}J
\quad\Longrightarrow\quad
\operatorname{ran}C^T=(\operatorname{im}J)^\perp.
\]

Consequently the neighbour pencil

\[
\begin{pmatrix}
A_{\rm loc}-zM_{\rm loc}&C^T\\
C&0
\end{pmatrix}

\]

has the same exact finite physical generalized spectrum as the assembled
Whitney pencil, by the already proved coefficientwise KKT argument.

- **DERIVED UNIFORMLY LOCAL SPECTRAL CONSTRAINT:** all gates pass at both
  levels with maximum degree at most three.
- **REFUTED:** any star graph disconnects, mixes simplices, exceeds the degree
  bound, or has the wrong exact kernel dimension.

Even a pass remains an algebraic descriptor system.  Its multiplier metric is
zero and redundant rows introduce multiplier gauge modes.  It is not a
positive-metric unitary tick, and no physical time claim is licensed.

Only the new targeted verifier will be run.
