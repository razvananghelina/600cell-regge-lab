# Prior-art gate: internal Regge-curvature response of the strong tangent modes

Date: 2026-08-17

## Exact object, carrier and hypotheses

This mission concerns the already-derived one-slab Lorentzian Regge solution,
not a continuum spacetime and not an arbitrary graph operator.

- The spatial boundary is the fixed 600-cell carrier with 120 vertices and 720
  edges on each boundary sheet.
- The four-dimensional slab has 2,280 edge variables: 720 old-boundary, 840
  internal and 720 new-boundary squared lengths.  Variations are logarithmic,
  so the coordinates are `delta log(q_e)` on the frozen Lorentzian branch.
- The slab has 6,240 triangle hinges.  The 2,400 boundary triangles are
  excluded from the primary observable because their `pi + sum(theta)` is a
  boundary/extrinsic-angle quantity.  The primary carrier is the 3,840
  internal triangles with deficit

  ```text
  epsilon_h = 2*pi + sum_{4-simplices s incident on h} theta_{s,h}.
  ```

- To place both Lorentzian causal types on a real branch, define the rowwise
  coordinate

  ```text
  kappa_h = -i epsilon_h  for a spacelike hinge,
  kappa_h =    epsilon_h  for a timelike hinge.
  ```

  This is an invertible fixed row phase.  Therefore it cannot create or remove
  a zero response.  Norms after this phase are **STRUCTURAL** coordinates, not
  a derived physical inner product on curvature.
- The observable is the Jacobian

  ```text
  D_kappa : C^2280 -> R^3840,
  (D_kappa)_{h,e} = partial kappa_h / partial log(q_e),
  ```

  evaluated at the committed homothetic dynamic solution.  It is assembled
  only from the already-audited local Lorentzian dihedral-angle formula and
  literal triangle--4-simplex incidence.
- The committed canonical boundary response supplies

  ```text
  Z : boundary phase C^1440 -> slab edge variations C^2280.
  ```

  The test operator is `F = D_kappa Z`.
- The frozen binary-tetrahedral action decomposes the calculation into seven
  deterministic minimal sectors.  From the target-blind tangent census, the
  strong expanding and contracting spaces have size `5d` in a sector of
  irrep dimension `d`, except that the trivial sector has size 4 rather than
  5 because its fifth reciprocal pair is the unresolved near-unit homogeneous
  pair.  Counting the regular-representation multiplicities gives exactly
  119 expanding and 119 contracting directions.
- Both independently frozen schedule parities and all four high-precision
  derivative variants remain mandatory.  No continuum harmonic, speed,
  Planck scale or desired curvature magnitude is loaded.

The hypotheses are therefore: the committed slab solution and branch are
correct; the canonical Legendre solve is the relevant one-step tangent; and
linear response of internal deficits is an admissible diagnostic of whether a
tangent direction changes Regge curvature.  None of these hypotheses says
that nonzero deficit response is sufficient to identify a physical graviton.

## Target-independent geometric census

Before looking at any strong tangent vector, a 100-decimal reconstruction gave
the same census for both schedule parities:

| hinge class | physical triangles | free `2T` orbits | base `kappa` |
|---|---:|---:|---:|
| boundary, spacelike | 2,400 | 100 | `+/- 0.000345939729225125...` |
| internal, spacelike | 2,400 | 100 | zero below `1e-40` |
| internal, timelike | 1,440 | 60 | `0.128388432031982800...` |

All 260 triangle orbits have size 24.  The minimum absolute triangle
area-square is `0.00020011381195962117733`, so there is no null-hinge branch
ambiguity.  After the fixed row phase, the maximum residual imaginary part of
the 6,240 base deficits is below `2.4e-96`.

This census is a **DERIVED COMPUTATIONAL CONTROL** of the carrier.  It is not a
result about the 119 tangent modes.

## Primary prior art

### KNOWN

1. In Regge calculus, curvature is concentrated on codimension-two hinges and
   is represented by deficit angles.  Hinge deficits, and deficits divided by
   appropriate dual areas when a scalar-curvature density is wanted, are
   standard Regge observables.  See McDonald and Miller,
   [*A Discrete Representation of Einstein's Geometric Theory of
   Gravitation*](https://arxiv.org/abs/0804.0279).
2. Around a flat Regge background, vertex displacements are exact gauge
   directions and gauge-invariant lattice gravitons can be organized as
   propagating curvature degrees of freedom.  See Hoehn,
   [*Canonical linearized Regge Calculus: counting lattice gravitons with
   Pachner moves*](https://arxiv.org/abs/1411.5672).
3. Exact discrete diffeomorphism symmetry is generically broken on curved
   Regge solutions.  The corresponding canonical relations become
   background-dependent pseudo-constraints rather than exact gauge
   constraints.  See Bahr and Dittrich,
   [*(Broken) Gauge Symmetries and Constraints in Regge
   Calculus*](https://arxiv.org/abs/0905.1670), and Dittrich and Hoehn,
   [*From covariant to canonical formulations of discrete
   gravity*](https://arxiv.org/abs/0912.1817).
4. Lorentzian angle conventions, including the complex continuation needed
   across causal sectors, require explicit branch control.  See Sorkin,
   [*Lorentzian angles and trigonometry including lightlike
   vectors*](https://arxiv.org/abs/1908.10022).

### CONTROL

- `D_kappa v = 0` is an exact, mechanically checkable statement that the
  linearized internal deficit vector is unchanged in direction `v`.
- `D_kappa v != 0` rules out membership in the kernel of this particular
  curvature observable.
- The result must be stable under both schedules, all derivative variants,
  the independent direct-eigenvector/ordered-Schur construction and the
  carried ball/roundoff uncertainty.

### OPEN

- No located primary source computes this exact 600-cell dust slab, its
  `2T`-resolved 119-dimensional strong tangent space or this response matrix.
  External novelty is **OPEN**, not claimed.
- On this curved background, zero response would support a pseudo-gauge
  interpretation but would not by itself prove a full gauge symmetry of the
  action.
- Nonzero response would show that the modes carry this Regge-curvature
  observable, but would not by itself distinguish physical gravitons from
  constraint-violating or scalar/vector mixtures.
- A Euclidean norm of the complete deficit vector is not a derived DeWitt or
  symplectic norm.  Only zero/nonzero, ranks, invariant subspaces and
  schedule/precision stability receive **DERIVED COMPUTATIONAL** status;
  relative gains remain **PATTERN** unless a physical norm is derived.
- Intrinsic curvature on the new three-dimensional boundary is a separate
  observable and remains a later gate.

## Framing attack

This test is necessary but not sufficient for the physical interpretation
suggested by the word “graviton.”  Deficit angles encode the full discrete
Riemann curvature at hinges, while Einstein evolution additionally involves
constraints and a Ricci/Weyl interpretation.  A large response must therefore
not be advertised as a gravitational wave.  Conversely, the curved-background
literature forbids the shortcut “small Legendre eigenvalue means gauge.”

The scientifically sharp question is only:

> Do the already-selected 119 expanding/contracting invariant tangent
> directions lie in the exact kernel of the internal deficit-angle Jacobian,
> or do they produce a resolved nonzero internal-curvature response?

That question is finite, falsifiable and independent of any desired continuum
dispersion relation.
