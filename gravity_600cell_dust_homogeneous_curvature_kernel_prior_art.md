# Prior-art gate: identity of the homogeneous curvature-kernel line

Date: 2026-08-17

## Why this mission is explicitly confirmatory

The preregistered internal-curvature census produced two facts before this
mission was formulated:

```text
rank(F) = 1439 out of 1440 boundary-phase directions,
rank(F restricted to either 119-dimensional strong branch) = 119.
```

The unique calibrated kernel lies in the trivial binary-tetrahedral sector.
The earlier target-blind tangent census also contains a unique non-strong real
reciprocal pair near `-1`, with moduli `0.9939037...` and `1.0061337...`.
Comparing those objects is therefore target-driven.  It can identify a line,
but it cannot be advertised as an unanticipated numerical hit.

## Exact object and carrier

For each frozen schedule parity, restrict the already-derived maps to the
60-dimensional trivial `2T` boundary-phase block:

```text
T : C^60 -> C^60                         canonical one-step tangent,
Y : C^60 -> C^65                         slab edge response,
F : C^60 -> C^160                        internal deficit response.
```

Coordinates are ordered as 30 old-boundary logarithmic squared-edge
positions followed by 30 conjugate momenta.  `Y` returns the 35 internal and
30 new-boundary logarithmic squared-edge variations.  `F = D_kappa Z` is the
incidence-derived real-branch deficit response certified in commit
`a3d9e7d`.

Let `K = ker(F)`, a calibrated one-dimensional right-singular subspace.  The
mission is to identify `K` against a finite, disclosed candidate list and to
test whether `K` is invariant under `T`.

## Candidate list fixed by existing geometry

The following candidates exhaust the comparisons in this mission:

1. the individual contracting near-`-1` eigenline of `T`;
2. the individual expanding near-`-1` eigenline of `T`;
3. their invariant two-plane;
4. the exact uniform boundary-position line `(1_30, 0_30)`;
5. the exact uniform boundary-momentum line `(0_30, 1_30)`;
6. their exact two-dimensional uniform phase plane;
7. the 30-dimensional pure-position subspace;
8. the 30-dimensional pure-momentum subspace;
9. after transport by `Y`, the five-dimensional canonical weak Schur-lift
   subspace selected by the five pole-edge orbit types;
10. after transport by `Y`, the independently derived five-dimensional
    geometric vertex-lapse subspace.

There are exactly 10 candidate comparisons per schedule and 20 total.  In
addition, one tangent-invariance test per schedule and one cross-schedule
kernel comparison are reported.  Combining candidates after reading their
distances or defining a new line from the observed kernel is forbidden.

Candidates 1--3 and the very decision to run this mission are
**post-observed**.  Candidates 4--10 are geometrically defined without using
the kernel vector, but they are still part of a disclosed 20-attempt census.

## Primary prior art

### KNOWN

- On flat Regge backgrounds, vertex displacements generate exact gauge
  directions and curvature variables separate the gauge and lattice-graviton
  sectors.  Hoehn,
  [*Canonical linearized Regge Calculus*](https://arxiv.org/abs/1411.5672).
- On curved Regge backgrounds, exact vertex-displacement symmetry is
  generically broken and canonical constraints become pseudo-constraints.
  Bahr and Dittrich,
  [*(Broken) Gauge Symmetries and Constraints in Regge
  Calculus*](https://arxiv.org/abs/0905.1670); Dittrich and Hoehn,
  [*From covariant to canonical formulations of discrete
  gravity*](https://arxiv.org/abs/0912.1817).
- Linearized Regge operators can have geometrically structured kernels tied
  to the discrete differential complex; a kernel must be interpreted through
  its carrier rather than by nullity alone.  Christiansen,
  [*On the linearization of Regge calculus*](https://arxiv.org/abs/1106.4266).

### CONTROL

- The right-singular kernel must have a resolved one-dimensional separation
  from the next singular value under all four derivative variants.
- Candidate subspaces must be fixed independently of the kernel vector.
- Equality means equality of subspaces within propagated finite-difference,
  Flint-ball, tangent-ball and binary-linear-algebra uncertainty; a small raw
  angle is not enough.
- A tangent eigenline identification requires both line equality and
  `T K = K` within the same calibrated gate.
- The two schedule kernels must be compared in the literal common
  old-boundary phase ordering, not after a fitted change of basis.

### OPEN

- No primary source found in the pre- or post-result searches contains this
  exact 600-cell dust kernel or these ten comparisons.  External novelty is
  **OPEN**.
- Membership in the geometric lapse subspace would identify the direction's
  kinematics, but on the curved background it would still not restore an
  exact diffeomorphism theorem for the action.
- Membership in the near-unit tangent plane would connect the count `119+1`
  to the curvature split, but would not show nonlinear integrability.
- Failure of all candidates would not remove the kernel; it would establish
  that its identity requires a new algebraic or geometric construction.

## Framing attack

Internal deficit angles are not a complete set of canonical constraints.  A
line in their kernel may arise from a discrete Bianchi dependence, a uniform
scale/momentum combination, a lapse-like deformation, or an accidental
linearized cancellation.  Therefore the strongest admissible conclusion is
“this line equals this preregistered subspace.”  The labels “time,” “gauge”
and “physical vacuum direction” remain **OPEN** until the relevant action
symmetry or nonlinear continuation is proved.
