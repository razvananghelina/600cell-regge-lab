# Dual incidence canonically resolves every constraint redundancy

Date: 2026-08-12

Preregistration commit: `c5f9bee`

Targeted verifier:
`reproducible/verify_whitney_dual_constraint_resolution.py`

Targeted result: **9/9 PASS**.  The verifier is registered exactly once.  No
spectrum, metric, independent-row choice or phenomenological target was
used.  The full suite was not run, by explicit user request.

## Headline

The redundant local Whitney copy constraints possess a complete canonical
resolution.  Relations among neighbour constraints are exactly boundaries of
dual two-cells, and relations among those relations are exactly boundaries of
dual three-cells.

> **DERIVED CANONICAL REDUCIBILITY RESOLUTION:** complete signed dual-cell
> incidence resolves every canonical neighbour-constraint relation at the
> base 600-cell and its first barycentric subdivision.

The construction removes the noncanonical spanning-tree choice that appeared
when one retained only an independent constraint basis.  It does not create
a physical gauge symmetry: the original copy constraints remain second
class.

The naive complete dual-cell maps also fail the preregistered bounded-degree
test across barycentric subdivision:

> **DERIVED NEGATIVE FOR THE NAIVE UNIFORMLY-BOUNDED HIERARCHY UNDER THIS
> BARYCENTRIC TRANSITION:** the largest dual-face boundary grows from 5 to 10
> incidences and the largest dual-volume boundary from 12 to 62.

This locality negative is specific to complete unsubdivided dual cells and
the tested barycentric transition.  The theory's later rank-edgewise tower is
a separate case and is not decided by this calculation.

## Canonical complex

For every global `p`-simplex `s`, its tetrahedron occurrences are the dual
vertices of the closed star of `s`.  A triangle containing `s` gives a dual
edge between its two parent tetrahedra.  When present, edges containing `s`
give dual faces, and a vertex gives one dual volume.

With `A_p` the unweighted quotient map from occurrences to their global
simplex, the verified sequence is

\[
 0\longrightarrow Z_{p,3}
 \mathrel{\mathop{\longrightarrow}^{R_{p,3}}}
 Z_{p,2}
 \mathrel{\mathop{\longrightarrow}^{R_{p,2}}}
 Z_{p,1}
 \mathrel{\mathop{\longrightarrow}^{C_p^T}}
 Z_{p,0}
 \mathrel{\mathop{\longrightarrow}^{A_p}}
 W_p\longrightarrow0.
\]

Every nonzero boundary coefficient is `+1` or `-1`.  Dual edges are oriented
by increasing tetrahedron index, dual face cycles by the least-index
convention, and the unique coherent dual-volume orientation has its
least-index face positive.  Changing an overall cell orientation only
multiplies a column by `-1` and cannot change exactness.

No metric, averaging weight, spanning tree or solver-selected null vector is
part of the definition.

## Complete dimensions and ranks

For a closed 3-complex with f-vector `(V,E,F,T)`, the preregistered flag
counts were

\[
\begin{array}{c|rrrr|r}
 &Z_0&Z_1&Z_2&Z_3&W_p\\
p=0&4T&3F&2E&V&V\\
p=1&6T&3F&E&0&E\\
p=2&4T&F&0&0&F\\
p=3&T&0&0&0&T.
\end{array}
\]

All are reproduced exactly.

### Base 600-cell

| degree | dimensions `(Z0,Z1,Z2,Z3,W)` | ranks `(A,C,R2,R3)` |
|---:|---:|---:|
| 0 | `(2400,3600,1440,120,120)` | `(120,2280,1320,120)` |
| 1 | `(3600,3600,720,0,720)` | `(720,2880,720,0)` |
| 2 | `(2400,1200,0,0,1200)` | `(1200,1200,0,0)` |
| 3 | `(600,0,0,0,600)` | `(600,0,0,0)` |

The second-stage rank is `1320+720=2040`, exactly the previously unexplained
row-redundancy count.  The 120 third-stage relations are the coherent
boundaries of the dual 3-cells at vertices.

### First barycentric subdivision

| degree | dimensions `(Z0,Z1,Z2,Z3,W)` | ranks `(A,C,R2,R3)` |
|---:|---:|---:|
| 0 | `(57600,86400,34080,2640,2640)` | `(2640,54960,31440,2640)` |
| 1 | `(86400,86400,17040,0,17040)` | `(17040,69360,17040,0)` |
| 2 | `(57600,28800,0,0,28800)` | `(28800,28800,0,0)` |
| 3 | `(14400,0,0,0,14400)` | `(14400,0,0,0)` |

Here the second-stage rank is

\[
 31440+17040=48480,
\]

exactly the refined redundancy count.

## Exactness certificate

The verifier checks all successive products as sparse integer matrices:

\[
 A_pC_p^T=0,
 \qquad C_p^TR_{p,2}=0,
 \qquad R_{p,2}R_{p,3}=0.
\]

Every product has exactly zero nonzero entries in all eight degree-level
blocks.

Degree-zero dual two-boundary ranks are computed blockwise modulo the two
primes `1000003` and `1000033`, and they agree.  The degree-one columns have
disjoint global-edge supports and are visibly independent; higher degrees
have no two-cells.  A nonzero minor modulo either prime is an exact lower
bound over the rationals, while nilpotency with the adjacent exact ranks gives
the matching upper bound.  Therefore these are exact rational ranks, not
numerical rank estimates.

Dimension equality then proves at every term:

\[
 \ker A_p=\operatorname{im}C_p^T,
\]

\[
 \ker C_p^T=\operatorname{im}R_{p,2},
\]

and

\[
 \ker R_{p,2}=\operatorname{im}R_{p,3}.
\]

Thus the sequence is a genuine resolution of conforming copies, not merely a
list of some visible cycle relations.

## Locality audit

The original neighbour boundary remains uniformly sparse:

- every constraint column has two occurrence endpoints;
- maximum occurrence-node row degrees remain `(3,2,1,0)`.

The higher complete dual-cell boundaries behave differently.

| map and degree | base maximum | barycentric maximum | factor |
|---|---:|---:|---:|
| quotient `A_0`, occurrences per vertex | 20 | 120 | 6 |
| quotient `A_1`, occurrences per edge | 5 | 10 | 2 |
| `R_0,2`, faces around an edge | 5 | 10 | 2 |
| `R_1,2`, faces around an edge | 5 | 10 | 2 |
| `R_0,3`, edges incident at a vertex | 12 | 62 | 5.1667 |

The refined cycle-length histograms are exact:

- `R_0,2`: 21,600 columns of length 4, 9,600 of length 6, and 2,880 of
  length 10;
- `R_1,2`: 10,800 columns of length 4, 4,800 of length 6, and 1,440 of
  length 10;
- `R_0,3`: 1,200 columns of length 8, 720 of length 12, 600 of length 14,
  and 120 of length 62.

The growth is not a floating support artefact: these are complete signed
incidence columns.

The negative does not yet apply to
`Esd_k(sd K)`.  Edgewise refinement was adopted precisely because iterated
barycentric refinement is not shape regular in the desired microscopic
sense.  The same resolution should next be tested on the preregistered
rank-edgewise controls before drawing an all-tower locality conclusion.

## Physical meaning and limitation

This result supplies a canonical hierarchy of multiplier relations.  In BRST
language it is the kinematic geometry needed for ghosts-for-ghosts when all
canonical neighbour rows are retained.

It does **not** change the exact bracket

\[
 CM_{\rm loc}^{-1}C^*,
\]

whose rank already proves that the independent copy constraints are second
class.  Nor does it supply an auxiliary Poisson structure, positive metric,
gauge-invariant Hamiltonian or physical state embedding.  The prior
first-class conversion result also remains intact: exact physical dressing
contains the global Gram inverse.

Therefore the honest gain is:

- the constraint redundancy is now completely derived and choice-free;
- the physical gauge and unitary-dynamics gates remain closed.

## Status ledger

- **DERIVED:** complete dual-cell relation maps with signed incidence only.
- **DERIVED:** exactness in every degree at both tested levels.
- **DERIVED:** ranks `2040` and `48480` explain all canonical-row
  redundancies.
- **DERIVED:** no spanning tree or independent-row basis is needed for the
  reducibility hierarchy.
- **DERIVED NEGATIVE:** complete unsubdivided dual cells do not retain
  bounded incidence under the tested barycentric refinement.
- **STRUCTURAL POSITIVE:** canonical kinematic input for a reducible BRST
  construction.
- **STRUCTURAL NEGATIVE:** this does not convert the original second-class
  constraints into physical gauge freedom.
- **OPEN:** locality on the actual rank-edgewise tower.
- **OPEN:** a geometry-selected auxiliary bracket, positive metric and local
  gauge-invariant Hamiltonian.
- **NOT CLAIMED:** a physical tick, time, causality, mass or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_dual_constraint_resolution.py
```

Expected result: `9/9`.
