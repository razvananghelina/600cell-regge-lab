# Preregistration: canonical dual-cell resolution of copy constraints

Date: 2026-08-12

This protocol is committed before constructing the new higher relation maps.
The row counts and first-stage ranks are already certified by
`verify_whitney_neighbour_constraints.py`; they are calibration inputs, not
new evidence.  No spectrum or phenomenological target is used.

## Question with complete hypotheses

Let `K` be either the 600-cell boundary triangulation or its first
barycentric subdivision.  It is a closed oriented combinatorial
three-manifold.  For each global `p`-simplex `s`, the duplicated Whitney
carrier contains one node `(T,s)` for every tetrahedron `T` containing `s`.

The already-derived neighbour constraint matrix `C_p` has one signed row for
every triangle `f` containing `s`; it equates the two copies in the two
tetrahedra adjacent across `f`.  Thus `C_p` is the signed incidence matrix of
the one-skeleton of the dual cell `D(s)`.

The new question is whether geometry canonically supplies the complete
relation hierarchy

\[
 0\longrightarrow Z_{p,3}
 \mathrel{\mathop{\longrightarrow}^{R_{p,3}}}
 Z_{p,2}
 \mathrel{\mathop{\longrightarrow}^{R_{p,2}}}
 Z_{p,1}
 \mathrel{\mathop{\longrightarrow}^{C_p^T}}
 Z_{p,0}
 \mathrel{\mathop{\longrightarrow}^{A_p}}
 W_p\longrightarrow0,
\]

where:

- `Z_p,0` consists of tetrahedron occurrences `(T,s)`;
- `Z_p,1` consists of flags `(s subset f)` with `f` a triangle;
- `Z_p,2` consists of flags `(s subset e)` with `e` an edge, when `p=0`,
  and one cell per global edge when `p=1`;
- `Z_0,3` consists of one cell per global vertex;
- layers above dual dimension `3-p` are zero;
- `A_p` sends every occurrence to its global simplex with its orientation
  sign, without averaging weights.

The proposed maps are the cellular boundary maps of `D(s)`.  Orient every
dual edge by increasing tetrahedron index.  For each dual two-cell, orient
its boundary cycle by the least-index convention; for each dual three-cell,
fix its unique coherent face orientation by setting its least-index face
coefficient positive.  Overall cell signs do not affect exactness or ranks.

This is a finite, label-deterministic construction.  There are no spanning
trees, metric weights, fitted coefficients or independent-row choices.

## Frozen exact gates

For every `p=0,1,2,3` at both levels:

1. reconstruct `C_p` and verify its existing exact kernel/rank certificate;
2. construct all dual two-cell boundary columns and require coefficient
   alphabet `{-1,+1}`;
3. construct all dual three-cell boundary columns and require the same
   alphabet;
4. verify exactly over the integers

   \[
   A_pC_p^T=0,\qquad C_p^TR_{p,2}=0,
   \qquad R_{p,2}R_{p,3}=0;
   \]

5. certify every rank blockwise over finite fields and use the displayed
   nilpotency bounds to turn matching lower bounds into exact rational ranks;
6. require exactness at every term, including

   \[
   \ker A_p=\operatorname{im}C_p^T.
   \]

The expected stage counts, derived from flags of the already-known
f-vectors `(V,E,F,T)`, are frozen as

\[
\begin{array}{c|rrrr}
 &Z_0&Z_1&Z_2&Z_3\\
p=0&4T&3F&2E&V\\
p=1&6T&3F&E&0\\
p=2&4T&F&0&0\\
p=3&T&0&0&0.
\end{array}
\]

Consequently the complete second-stage counts should be `3E` and the
third-stage count `V`.  Their ranks are not assumed; they are part of the
test.

## Frozen locality audit

Record for every map:

- maximum nonzeros per column and row;
- the complete histogram when practical;
- the ratio of each maximum between the base and barycentric level.

Decision labels:

- exactness with purely signed incidences is a **DERIVED CANONICAL
  REDUCIBILITY RESOLUTION**;
- if every maximum incidence degree stays unchanged, report only a
  **STRUCTURAL BOUNDED-LOCALITY POSITIVE ON TWO LEVELS**, not an all-level
  theorem;
- if any maximum grows, report a **DERIVED NEGATIVE FOR THE NAIVE
  UNIFORMLY-BOUNDED HIERARCHY**.  A canonical subdivision of the offending
  dual cells remains a distinct future construction.

No map may be discarded because it has high degree; complete incidence is
the preregistered canonical rule.

## Physical framing attack

Even a successful resolution does not turn the copy constraints into a
physical gauge symmetry.  Their bracket on the original field carrier is
already proved second class.  The new maps would resolve the *redundancy of
the full canonical row set* and remove the noncanonical independent-row
choice.  They would provide kinematic input for a reducible BRST extension,
not select:

- an auxiliary symplectic bracket;
- a positive metric;
- a gauge-invariant Hamiltonian;
- a physical embedding or a time scale.

The minimal first-class conversion already showed that exact physical
dressing introduces a global Gram inverse.  A local relation hierarchy does
not refute that result.

## Acceptance and kill boundaries

Acceptance for this sub-route requires exactness of the complete hierarchy
with maps selected only by dual incidence.  This advances the kinematic
geometry but not the physical tick.

If the incidence maps fail exactness, the proposed canonical reducible
hierarchy is killed.  If exactness succeeds but locality degree grows, the
naive dual-cell hierarchy is retained as a finite exact resolution but is
killed as a uniformly bounded microscopic rule.

## Outputs and exclusions

The registered verifier will write
`reproducible/whitney_dual_constraint_resolution.json`.  The result note will
include the protocol commit, full counts, ranks, incidence degrees and the
physical limitation.

Excluded:

- no spectrum or target comparison;
- no metric or BRST Hamiltonian;
- no claim that multiplier redundancy is physical gauge freedom;
- no time, causality, mass or Planck units;
- no full-suite run, by explicit user request.
