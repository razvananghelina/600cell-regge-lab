# Protocol: exact feasibility census for the refined canonical map

Date: 2026-08-20

Prior-art gate: commit `883b4e7`.

This protocol is frozen before constructing any refined action Hessian and
before inspecting any refined spectrum.

## 1. Frozen inputs

The verifier must hash and reconstruct from:

```text
commons/cell600.py
reproducible/verify_gravity_600cell_projected_rank_edgewise_carrier.py
reproducible/gravity_600cell_projected_rank_edgewise_carrier.json
reproducible/verify_gravity_600cell_projected_rank_edgewise_balanced_slab.py
reproducible/gravity_600cell_projected_rank_edgewise_balanced_slab.json
```

The only carriers are

```text
K0 = P(sd K_600),
K1 = P(Esd_2(sd K_600)).
```

The base colour is face dimension.  The fine colour of an edgewise vertex
key `(v,w)` is `c(v)+c(w) mod 4`.  No colour order is privileged.

## 2. Exact staircase convention

For each of all `4!` colour orders, sort the vertices of every spatial
tetrahedron as `v0<v1<v2<v3`.  The four standard product pentachora are

```text
{(v0,0),...,(vk,0),(vk,1),...,(v3,1)},  k=0,1,2,3.
```

Use literal labelled layer vertices.  An internal edge is any slab edge not
contained entirely in the old or new spatial boundary.

The verifier must derive rather than assume the identities

```text
N4_slab       = 4*T,
E_boundary    = 2*E,
E_internal    = V+E,
E_total       = V+3*E,
dim Legendre  = E_internal+E = V+2*E,
dim phase map = 2*E.
```

It must check every identity independently against explicit simplex and edge
sets for every schedule at both levels.

## 3. Schedule census

For every order record the internal-edge digest.  Report:

1. the number of distinct internal-edge sets;
2. the intersection and union sizes of cross-layer nonvertical edges;
3. the six spatial colour-pair edge populations;
4. for all `binomial(24,2)=276` schedule pairs, the weighted inversion
   distance: the number of spatial edges whose chosen cross diagonal changes;
5. the complete distance histogram, minimum positive distance and maximum.

As an independent formula control, the weighted inversion distance must equal
the sum of the colour-pair populations over pairs whose relative order is
reversed.  Direct cross-edge set comparison and the formula must agree for all
276 pairs.

## 4. Size and sparsity ledger

For each level print exact decimal and binary storage sizes for dense
`float64` arrays representing:

```text
the complete slab Hessian,          dimension E_total;
the pre-Legendre Jacobian,          dimension V+2E;
one boundary canonical phase map,  dimension 2E.
```

These sizes are facts, not machine-independent impossibility theorems.

A pentachoron has ten edges and therefore contributes at most

```text
10*11/2 = 55
```

upper-triangular Hessian positions.  Report the exact local-incidence upper
bound `55*N4_slab`.  Do not materialize a fine dense matrix or claim that this
upper bound equals the number of distinct sparse entries.

## 5. Controls and exclusions

Mandatory controls:

- both frozen spatial f-vectors, Euler characteristics and face incidence;
- proper four-colouring and one occurrence of each colour per tetrahedron;
- exactly 24 orders and the frozen pentachoron counts;
- exactly two spatial boundary components and no side boundary, inherited
  from the upstream slab certification;
- formula/direct equality for all edge and distance counts;
- byte-identical rerun of the final JSON artifact.

Forbidden in this mission:

- any Regge derivative or Hessian value;
- any eigenvalue, mode, continuum target, speed or Planck quantity;
- selecting or averaging a schedule;
- constructing a coarse-to-fine transport;
- interpreting a memory size as a physical no-go.

## 6. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_MAP_FEASIBILITY_CONTROL_FAILED` if provenance, topology,
   colouring, staircase or direct/formula controls fail.
2. `REFINED_MAP_SINGLE_TEMPORAL_CARRIER` only if the declared geometry yields
   one distinct internal-edge set after fixed time orientation.
3. `REFINED_MAP_SCHEDULE_ELIMINATION_REQUIRED` if multiple legitimate
   schedule edge sets remain.

Outcome 3 licenses only an all-schedule effective-boundary covariance test on
`K0`, using sparse or matrix-free elimination.  It does not license choosing
one order on `K1`.  Spatial refinement comparison remains blocked until both
temporal schedule dependence and an inter-level phase-space transport are
settled independently.

## 7. Deliverables

- registered verifier
  `reproducible/verify_gravity_600cell_refined_canonical_map_feasibility.py`;
- deterministic JSON artifact of the same stem;
- result note with the exact counts and scoped verdict;
- targeted verifier plus static registry/coverage checks only; no full suite.
