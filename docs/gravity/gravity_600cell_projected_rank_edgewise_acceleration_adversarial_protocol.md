# Adversarial protocol: orbit-compressed canonical-carrier acceleration

Date: 2026-08-19

Primary blind commit: `9469e33`  
Primary comparison commit: `f65dfc8`

## Independence target

The primary calculation applies one vectorized action implementation to every
tetrahedron and obtains logarithmic derivatives through a mixture of complex
and real centred differences.

The audit must not reuse that decisive numerical path.  It will:

1. reconstruct the 600 source tetrahedra with maximal-clique enumeration;
2. rebuild the projected barycentric carrier and its direct rank-selected
   eight-child refinement;
3. canonically classify unique spatial edges, faces and tetrahedra by their
   intrinsic distance data up to vertex permutation;
4. evaluate one representative per intrinsic shape class and multiply by the
   independently counted class multiplicity;
5. use a real five-point logarithmic derivative stencil, different step
   sizes, different `eta` values and different coefficient sentinels from the
   primary calculation.

Thus the audit replaces more than 100,000 simultaneous metric inversions by a
small weighted shape census and replaces the derivative mechanism.  It uses
the same declared Regge action; changing the physical action would not be an
independent verification of the same claim.

## Frozen provenance

```text
primary blind artifact
  reproducible/gravity_600cell_projected_rank_edgewise_acceleration_blind.json
  2059620f22cfbd8eac8abe6f2c7536924128d37f47a430bf773e34a9aead93a2

primary comparison artifact
  reproducible/gravity_600cell_projected_rank_edgewise_acceleration_comparison.json
  132a81fe03ee67dbe95b91a68910f9212db88f8c6104b23f7d3f3f422939f5a4

canonical carrier artifact
  reproducible/gravity_600cell_projected_rank_edgewise_carrier.json
  b57955b85a972df00b5673ddf7ee295757848f5afb43314857cf3de2dc85ac84
```

## Frozen intrinsic classification

- Edge signature: squared chord length rounded to 11 decimal places.
- Face signature: the sorted three squared chord lengths, each rounded to 11
  decimal places.
- Tetrahedron signature: lexicographically least six-entry squared-distance
  tuple over all 24 vertex permutations, rounded to 11 decimal places.

For every resulting class, the maximum deviation from its representative's
unrounded signature must be `<2e-9`.  Edge, face and tetrahedron
multiplicities must sum to the complete frozen f-vector at each level.

Grouping is used only after this within-class residual gate.  No number of
classes or multiplicity is assumed in advance.

## Frozen reduced action

For one representative of each tetrahedron class, compute directly:

- its spatial Gram matrix and intrinsic spatial dihedral angles;
- its Lorentzian frustum Gram matrix;
- all six lateral angles and all eight boundary angles;
- local lateral-hinge and boundary-face contributions.

The global `2*pi` edge and `pi` boundary-face terms are summed from the edge
and face class censuses.  Every representative local term is multiplied by
its tetrahedron-class multiplicity.  Dust mass and volume-radius
normalization are selected exactly as in the primary protocol.

At one static and two non-static held-out states, the compressed action must
agree with a complete uncompressed evaluation to relative error `<2e-8`.
These held-out comparisons validate compression, not the final coefficient.

## Different derivative and extrapolation controls

Use the real five-point stencil

```text
(f(-2h)-8*f(-h)+8*f(h)-f(2h))/(12h)
```

in logarithmic coordinates, with

```text
seam eta values       (0.05,0.025,0.0125,0.00625)
lapse eta values      (0.05,0.025,0.0125)
coefficient sentinels (0,-0.75,-1.5,-2.25)
seam log steps        4e-4 and 2e-4
lapse log steps       4e-3 and 2e-3.
```

Use the same algebraically required affine seam root, quadratic lapse root
after removal of the static root, and fourth-order Richardson hierarchy.  The
audit coefficient routes must agree with each other within `2e-5`, and each
audit coefficient must agree with its frozen primary coefficient within
`2e-5`.

## Positive and negative controls

Positive control: the orbit-compressed implementation on the regular
600-cell must reproduce the already certified radius coefficient

```text
-0.5394897340206755
```

within `2e-5`.

Negative control: change the multiplicity of the largest tetrahedron class by
one while leaving every global class census fixed.  At the held-out
non-static state this deliberately inconsistent action must differ from the
correct compressed action by more than `1e-8` relative.  The purpose is to
show that the compression check can detect a wrong weighted census.

## Outcome hierarchy

1. `ADVERSARIAL_CANONICAL_ACCELERATION_DISAGREEMENT` if provenance,
   intrinsic-class, compression, positive/negative-control or coefficient
   agreement fails.
2. `ADVERSARIAL_CANONICAL_ACCELERATION_CORROBORATED` only if every gate
   passes.

A corroborated outcome remains **DERIVED NUMERICAL CONTROL**, not a new law
of gravity.  One canonical refinement does not prove asymptotic order.
