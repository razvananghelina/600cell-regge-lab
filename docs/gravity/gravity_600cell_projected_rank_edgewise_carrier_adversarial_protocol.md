# Adversarial protocol: direct-split projected rank-edgewise carrier

Date: 2026-08-19

Independence gate commit: `88b12a0`

Status: frozen before constructing or hashing the direct-split carrier.

## 1. Frozen input

Require commit `cb4fc24`'s artifact

```text
reproducible/gravity_600cell_projected_rank_edgewise_carrier.json
SHA-256 b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84
```

with outcome `PROJECTED_RANK_EDGEWISE_CARRIER_DERIVED` and `16/16` passing
checks.  Require the source `commons/cell600.py` hash already fixed by the
primary protocol.

## 2. Independent source and base construction

Use `networkx.find_cliques` on the 600-cell graph.  Require exactly 600
maximal cliques, all of size four.  Construct all nonempty cells and complete
containment flags independently.  Reproduce

```text
f(K)=(120,720,1200,600),
f(P(sd K))=(2640,17040,28800,14400).
```

The vertex position for a cell remains the normalized mean of its original
600-cell vertices because that is the declared mathematical object, not an
implementation choice.

## 3. Direct eight-child split

For every rank-ordered chamber `(v0,v1,v2,v3)`, create global old-vertex and
edge-midpoint keys and use exactly

```text
(v0,m01,m02,m03)       (v1,m01,m12,m13)
(v2,m02,m12,m23)       (v3,m03,m13,m23)

(m01,m02,m03,m13)      (m01,m02,m12,m13)
(m02,m03,m13,m23)      (m02,m12,m13,m23).
```

This direct list is fixed before execution.  It must yield globally

```text
f=(19680,134880,230400,115200),
```

with no duplicate top cell, face incidence exactly two and Euler
characteristic zero.

Store a SHA-256 digest of the lexicographically sorted integer top-simplex
array and a separate digest of all projected vertex coordinates serialized
as little-endian binary64.  These are audit provenance, not preregistered
target hashes.

## 4. Local canonicity and negative control

On an abstract parent tetrahedron:

1. union the direct eight-child list over all 24 rank-labelled barycentric
   chambers;
2. require 192 distinct children and invariance under all 24 parent
   permutations;
3. apply the same central diagonal `m02--m13` directly to the unranked parent;
4. among the 12 even permutations require at least one preservation and at
   least one failure, with fewer than 12 preservations.

The last case must fail full `A4` invariance and is the deliberately wrong
control.  A control failure invalidates the audit rather than weakening its
conclusion.

## 5. Independent geometry

For every actual tetrahedron compute all six squared edge lengths and obtain
the volume only from the `5 x 5` Cayley--Menger determinant

```text
288 V^2 = det(CM).
```

Use it to reconstruct the same quality, chord, sag and total-volume fields as
the primary artifact.  First require on a regular tetrahedron that
Cayley--Menger and a separate Gram formula agree relatively below `1e-13`.

For both levels, every scalar field in

```text
quality_minimum, quality_median, quality_maximum,
volume_minimum, volume_maximum, volume_total,
volume_target_2pi2, volume_absolute_error,
maximum_chord_length, maximum_centroid_radial_sag
```

must agree with the frozen primary value within

```text
5e-11 * max(1,abs(primary value)).
```

The audit must independently satisfy the primary shape/chord/sag inequalities;
agreement with a favorable stored scalar is not enough if an inequality fails.

## 6. Outcome hierarchy

- `ADVERSARIAL_PROJECTED_RANK_EDGEWISE_CARRIER_CORROBORATED` if frozen
  provenance, source/base reconstruction, direct global topology, local
  positive/negative controls, Cayley--Menger calibration and every geometry
  comparison pass;
- `ADVERSARIAL_PROJECTED_RANK_EDGEWISE_CARRIER_DISAGREEMENT_OPEN` if controls
  pass but the actual direct carrier disagrees;
- `ADVERSARIAL_PROJECTED_RANK_EDGEWISE_CARRIER_CONTROL_FAILED` otherwise.

Write

```text
reproducible/gravity_600cell_projected_rank_edgewise_carrier_adversarial.json
```

and require two byte-identical complete executions before consolidation.

## 7. Scope

Corroboration is **DERIVED COMPUTATIONAL** for the finite carrier.  Radial
selection of the round background remains **STRUCTURAL**.  All-level nonlinear
projection, dynamics, constraints, dust and continuum physics remain
**OPEN**.
