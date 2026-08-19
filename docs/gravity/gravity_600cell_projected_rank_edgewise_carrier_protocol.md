# Protocol: canonical projected rank-edgewise 600-cell carrier

Date: 2026-08-19

Prior-art gate commit: `b361eca`

Status: frozen before constructing the full projected `Esd_2(sd K)` carrier
or evaluating any new shape, sag or volume number.  The base-control wording
was corrected before implementation because the old JSON stores no coordinate
array; see Section 3.

## 1. Frozen inputs

Require the following SHA-256 values:

```text
commons/cell600.py
ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f

reproducible/verify_whitney_rank_edgewise_refinement.py
371e28210fcf54f87acba114d26e6ffa8e72246842ef75000b832b4a6847e5dc

reproducible/whitney_rank_edgewise_refinement.json
af0c615a104ee7e22e0f003bde249d90a1b661b5e3cabec242a5763cc44aa77f

reproducible/verify_smooth_hopf_refinement_blind.py
2c0d0cb2ee1a9b1b6e0d3df7c35031f07e49a16850ecd02c8536b5b371f92b8a

reproducible/smooth_hopf_refinement_blind.json
7258a4755ac32af9d32d2415c09bf65f7b1c4a064475a0e2da304f43eb362ba8
```

The earlier artifacts are controls only.  Reconstruct the new full carrier;
do not import a stored mesh.

## 2. Source complex and group control

Reconstruct the 120 unit quaternions and all 600 four-cliques of the
600-cell.  Require

```text
f(K)=(120,720,1200,600),
degree=12,
each triangle has two incident tetrahedra,
chi(K)=0.
```

Build the complete quaternion multiplication table by matching products back
to the 120 source vertices.  Require closure and maximum matching residual
below `2e-8`.

For one tetrahedron, enumerate the rotational actions

```text
q -> a q b^-1,       a,b in 2I,
```

that preserve its vertex set.  Their distinct induced permutations must be
exactly the 12 even permutations.  Adding quaternion conjugation must give
all 24 permutations.

Let the three possible central-octahedron diagonals be the three partitions
of four vertices into two unordered pairs.  Require that the rotational
stabilizer has one orbit of size three and fixes none.  This is the explicit
full-carrier control for the red-refinement no-go; it is not inferred merely
from an abstract tetrahedron.

## 3. Choice-free base `K_0=P(sd K)`

Index every nonempty face of `K`.  Place its new vertex at the normalized
Euclidean mean of its original 600-cell vertices.  A top simplex is every
complete flag

```text
vertex subset edge subset triangle subset tetrahedron.
```

The expected f-vector, fixed combinatorially, is

```text
f(K_0)=(2640,17040,28800,14400).
```

Require no duplicate top simplex, incidence two for every triangle and Euler
characteristic zero.  Reproduce the old artifact's stored dimensions and
maximum chord length `0.385707678423` within `2e-12`.  The coordinates
themselves are reconstructed by the formula frozen above and in the
source-hashed old verifier; the old JSON does not contain a coordinate array,
so a claimed direct array comparison would be impossible.  This correction
was made before implementing or running the new verifier and changes no new
carrier, shape threshold or outcome rule.

## 4. Rank-selected `Esd_2` and projection

Generate the local Edelsbrunner--Grayson color schemes independently at
`k=2`.  On one rank-ordered chamber require exactly eight tetrahedra, ten
vertices, equal affine volume `1/8` and no duplicate top cell.  The union over
all 24 complete flags of an abstract parent tetrahedron must contain exactly
192 tetrahedra and be invariant under all 24 parent-vertex permutations.

For every ordered chamber of `K_0`, substitute its four vertices in rank
order `(0,1,2,3)`.  Merge a weight-two endpoint with the corresponding old
vertex and merge a `(1,1)` weight with the midpoint of the same global edge.
Only after global merging, radially normalize every point.

The expected global f-vector is fixed before execution:

```text
f(K_1)=(19680,134880,230400,115200).
```

Here `19680=2640+17040`, `115200=8*14400`, closed face incidence gives
`230400=2*115200`, and Euler zero fixes `134880`.

Require all four counts independently, no duplicate tetrahedron, incidence
two for every triangle, Euler zero and unit vertex norms within `2e-12`.

## 5. Equivariance gate

For every one of the 120 left multiplications, 120 right multiplications and
quaternion conjugation:

1. map every original face to its indexed image;
2. require the complete `K_0` top-simplex set to be preserved;
3. compare the transformed projected face barycentre with the projected image
   face barycentre.

The maximum coordinate residual must be below `2e-8`.  Together with the
exact local `S4` edgewise test and rank preservation, this certifies the
global projected construction; an exhaustive `14400 x 115200` materialized
permutation table is unnecessary and would not add a new logical condition.

## 6. Frozen finite-geometry gates

For every straight-chord tetrahedron compute its three-volume `V`, all six
edge lengths and mean-ratio quality

```text
q = 12*(3V)^(2/3) / sum_edges(length^2).
```

For `K_0` and `K_1` record:

- minimum, median and maximum `q`;
- minimum and maximum tetrahedral volume;
- maximum chord length `h`;
- maximum element-centroid radial sag
  `s=max_t (1-||mean(vertices(t))||)`;
- total chordal volume and its absolute distance from `2*pi^2`.

The carrier passes the finite-geometry gate only if:

```text
all V > 0,
q_min(K_0) > 0.25,
q_min(K_1) > 0.25,
q_min(K_1) >= 0.5*q_min(K_0),
h(K_1) < h(K_0),
s(K_1) < s(K_0).
```

The `0.25` and factor-two gates are conservative finite-element hygiene
bounds fixed without seeing `K_1`.  Volume-error monotonicity is reported and
mechanically labelled `IMPROVES` or `DOES_NOT_IMPROVE`, but it is not allowed
to override topology, symmetry or shape.  One volume scalar is not a proof of
metric convergence.

## 7. Outcome hierarchy

- `PROJECTED_RANK_EDGEWISE_CARRIER_DERIVED` if every provenance, group,
  topology, local-choice, equivariance and finite-geometry gate passes;
- `PROJECTED_RANK_EDGEWISE_CANONICITY_FAILED` if a group, local-choice,
  topology or equivariance gate fails;
- `PROJECTED_RANK_EDGEWISE_FINITE_GEOMETRY_FAILED` if canonicity passes but a
  volume/quality/chord/sag gate fails;
- `PROJECTED_RANK_EDGEWISE_CONTROL_FAILED` for a frozen input or source
  calibration failure.

Passing establishes a choice-free first projected edgewise refinement.  It
does not establish all-level shape regularity after nonlinear projection,
Regge convergence or physical dynamics.

## 8. Artifact and falsification

Write every count and diagnostic to

```text
reproducible/gravity_600cell_projected_rank_edgewise_carrier.json
```

Two complete executions must be byte-identical.  Record the no-red orbit,
the complete local permutation ledger and both level diagnostics, not only
the final label.

Before consolidation, a second registered verifier must avoid the primary
color-scheme construction at the decisive first refinement.  It will rebuild
the `k=2` carrier directly from old vertices, global edge midpoints and an
independent local Freudenthal/pulling description, and compare canonical
simplex and geometry hashes.  It must include a deliberately non-equivariant
red-diagonal control.  Disagreement leaves the result **OPEN** under Rule 4.

## 9. Interpretation boundary

- **DERIVED** may label exact topology, stabilizer and equivariance facts.
- **STRUCTURAL** labels the use of the unit round sphere and radial map.
- **PATTERN** labels any one-step approximation trend.
- **OPEN** remains the infinite projected tower, a Lorentzian action,
  constraint reduction, graviton dispersion, local dust, speed and scales.

No continuum dispersion, polarization count, speed, Planck quantity or
particle target is loaded.
