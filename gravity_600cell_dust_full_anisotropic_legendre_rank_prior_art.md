# Prior-art gate: full anisotropic canonical rank of one 600-cell dust slab

Date: 2026-08-17

Status: written before evaluating any entry or singular value of the full
`1560 x 1560` canonical Jacobian.

## Exact object and complete hypotheses

Use the already accepted first non-static homogeneous dust slab, whose frozen
artifact has SHA-256

```text
4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9.
```

For each of the independently derived `even` and `odd` five-stage schedules,
retain every individual Regge edge rather than tying the order-24 orbits:

```text
old boundary spacelike edges       720
internal diagonals                 720
internal timelike poles            120
new boundary spacelike edges       720
complete squared-edge variables   2280.
```

The slab consists of 2400 Lorentzian four-simplices.  The old and new
boundaries are complete regular-600-cell edge sets.  Positive magnitudes are
used for the 120 poles while their actual squared lengths carry the inherited
minus sign.

Use logarithmic coordinates for all positive magnitudes.  Split the complete
action Hessian according to

```text
O = old boundary,  dim 720,
X = internal,      dim 840,
N = new boundary,  dim 720.
```

With the same action and momentum signs already certified in the reduced
calculation, the full pre-Legendre Jacobian is

```text
J_full = [[ K_XX,  K_XN],
          [-K_OX, -K_ON]],             shape 1560 x 1560.
```

It differentiates the equations that reconstruct `(X,N)` from fixed old
geometry and fixed old canonical momentum.  This mission asks only for its
rank, calibrated small-singular-value structure and boundary support of any
nullspace.  It does not yet construct or diagonalize a gravitational-wave
propagator.

Dust is a fixed conserved total mass divided uniformly among the 120 vertex
world-lines:

```text
S_dust = -(8*pi*M/120) sum_poles sqrt(rho_p).
```

At the homogeneous point this is exactly the already certified reduced dust
term.  The equal per-vertex localization is a stated physical hypothesis; a
different inhomogeneous matter discretization is not covered.

## Geometry-only reduction known before the Hessian

The ordered schedule stabilizer has order 24.  Direct enumeration before any
new Hessian calculation gives element orders

```text
1^1, 2^1, 3^8, 4^6, 6^8
```

and seven conjugacy classes of sizes

```text
1, 1, 4, 4, 4, 4, 6.
```

This identifies the stabilizer as the binary tetrahedral group `2T`.  Every
one of the 95 edge orbits is free and has size 24.  Consequently the full
2280-coordinate representation is 95 copies of the regular representation,
and `J_full` splits into the seven complex irreducible sectors of dimensions

```text
1,1,1,2,2,2,3.
```

The corresponding reduced canonical block sizes are

```text
65,65,65,130,130,130,195,
```

with each block spectrum repeated by its irrep dimension in the full matrix.
This finite symmetry decomposition is geometry-derived and target-free.  It
does not assert that the complete slab has the full `H4` symmetry; the chosen
schedule retains only `2T`.

## Primary prior art

Barrett, Galassi, Miller, Sorkin, Tuckey and Williams give a local implicit
Regge evolution scheme and illustrate it on a 600-cell Friedmann cosmology:
[*A Parallelizable Implicit Evolution Scheme for Regge
Calculus*](https://arxiv.org/abs/gr-qc/9411008).

De Felice and Fabri enlarge the set of variables in a Sorkin evolution of a
dust-filled 600-cell and study its causal endpoint:
[*Singularities of the closed RW metric in Regge Calculus: a generalized
evolution of the 600-cell*](https://arxiv.org/abs/gr-qc/0106077).

Dittrich and Hoehn derive pre/post Legendre evolution directly from the
simplicial action:
[*Canonical simplicial gravity*](https://arxiv.org/abs/1108.1974).
Their covariant-to-canonical analysis shows that flat-background vertex
displacements give exact constraints, whereas nonlinear curved backgrounds
generically produce background-dependent pseudo-constraints:
[*From covariant to canonical formulations of discrete
gravity*](https://arxiv.org/abs/0912.1817).

Hoehn's linearized Pachner-move analysis counts lapse/shift and lattice
graviton data on flat backgrounds only after the constraint generators are
identified:
[*Canonical linearized Regge Calculus: counting lattice gravitons with
Pachner moves*](https://arxiv.org/abs/1411.5672).
Bahr and Dittrich likewise stress that curvature can break exact discrete
diffeomorphism symmetry and replace constraints by pseudo-constraints:
[*\(Broken\) Gauge Symmetries and Constraints in Regge
Calculus*](https://arxiv.org/abs/0905.1670).

## KNOWN / CONTROL / OPEN

### KNOWN

- Action-generated pre/post Legendre maps and constrained canonical Regge
  evolution are established theory.
- Exact vertex-displacement gauge directions are expected on suitable flat
  backgrounds; they cannot be assumed on the present curved dust background.
- Sorkin-type 600-cell evolution and enlarged-variable dust calculations are
  prior art.
- The repository's order-24 invariant quotient has a calibrated full-rank
  `65 x 65` Jacobian at this same dynamic state.  That proves only the trivial
  `2T` sector and cannot rule out nulls in the six nontrivial sectors.

### CONTROL

- The full gradient must be constant on every free 24-element orbit and its
  orbit restriction must reproduce the certified 95-coordinate gradient.
- Restricting the assembled full Hessian to orbit-constant variations must
  reproduce the committed trivial-sector canonical matrix and tangent map.
- The Hessian must be reciprocal, the complete internal equations must vanish
  edge by edge, and all local simplex evaluations must stay on the certified
  Lorentzian branch.
- Local angle derivatives must agree across independent high-precision step
  pairs, and deterministic full-space directional derivatives must reproduce
  the assembled Hessian.
- The seven symmetry blocks must exhaust dimension 1560, reproduce the full
  singular multiset with the required irrep repetitions, and show negligible
  off-block leakage relative to their numerical calibration.

### OPEN

- Whether `J_full` is regular outside the invariant quotient.
- Whether small nonzero modes form calibrated pseudo-constraints or exact
  nullspaces.
- Whether any nullspace reaches the 720-dimensional new-boundary geometry.
- The full canonical tangent map, its physical constraint quotient and its
  stability.
- Tensor-mode identification, a discrete dispersion relation, refinement
  stability, a causal cone and a limiting speed.
- External novelty of this exact `2T`-resolved rank census.  A primary search
  located the general frameworks above but no identical matrix calculation;
  absence from a search is not proof of novelty.

## Framing attack

A full-rank result would not prove that all 720 boundary shapes are physical.
On a curved finite Regge background it may instead show that discretization
has lifted continuum gauge into small pseudo-constraints.  Calling every
resolved direction a graviton would therefore be wrong.

A null result would also not automatically identify gauge: a null vector must
be matched to an explicit vertex-displacement generator and its projection on
the final boundary must be measured.  Symmetry, an accidental degeneracy and
gauge are distinct hypotheses.

Finally, the selected staircase itself breaks `H4` to `2T`.  A spectrum from
one schedule is not a continuum tensor spectrum.  This rank census is only
the necessary gate deciding whether the complete propagator can be built at
all.
