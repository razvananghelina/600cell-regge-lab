# Prior-art gate: three-slice Jacobi operator of the dust 600-cell

Date: 2026-08-18

## Exact object, carrier and hypotheses

Use only the two already accepted, canonically sewn, fixed-total-mass Regge
dust slabs and their committed complete tangent balls

```text
T_1 : (delta q_0, delta p_0) -> (delta q_1, delta p_1),
T_2 : (delta q_1, delta p_1) -> (delta q_2, delta p_2).
```

Here every `q_n` consists of the `720` logarithmic signed-squared spatial
edge variables on boundary slice `n`; every `p_n` is its action-derived
canonical momentum.  The full phase dimension is `1440`.  The carrier is
resolved in the seven already derived binary-tetrahedral minimal sectors,
with position dimensions

```text
n = 90, 60, 60, 60, 30, 30, 30
```

and representation weights `3,2,2,2,1,1,1`.

For a tangent block

```text
T_i = [ A_i  B_i ]
      [ C_i  D_i ],
```

the proposed object is the linearized, internal-edge-eliminated seam equation
on three consecutive positions,

```text
K_- delta q_0 + K_0 delta q_1 + K_+ delta q_2 = 0.
```

When the discrete Legendre twist blocks `B_i` are invertible, reconstruct the
quadratic Hamilton principal functions from the maps:

```text
S_i,01 = -B_i^-1,
S_i,00 =  B_i^-1 A_i,
S_i,10 =  C_i - D_i B_i^-1 A_i,
S_i,11 =  D_i B_i^-1,

K_- = S_1,10,
K_0 = S_1,11 + S_2,00,
K_+ = S_2,01.
```

This gives the normalized recurrence

```text
delta q_2 = P delta q_1 + Q delta q_0,
P = -K_+^-1 K_0,
Q = -K_+^-1 K_-.
```

`Jacobi operator` here means the Hessian/linearized discrete Euler--Lagrange
operator along the solved trajectory.  It does not refer to the Jacobi
identity and does not assume a continuum wave equation.

All four committed derivative variants and both independently derived
schedule parities are mandatory.  The calculation must propagate the stored
Flint ball radii and the half-ULP enclosure of binary serialization.  It may
not fit coefficients or load a spatial spectrum.

## Primary prior art

- [Marsden--West](https://doi.org/10.1017/S096249290100006X) derive discrete
  Euler--Lagrange evolution, discrete Legendre transforms and symplectic maps
  from a discrete action.  Linearizing the discrete Euler--Lagrange equation
  gives the general three-point Jacobi recurrence used here.
- [Dittrich--Hoehn](https://arxiv.org/abs/1108.1974) formulate Regge evolution
  as action-generated composable canonical transformations, including
  changing hypersurfaces and constraints.
- [Dittrich--Hoehn](https://arxiv.org/abs/0912.1817) derive linearized
  canonical Regge dynamics and explain its background-dependent
  pseudo-constraints.
- [Hoehn](https://arxiv.org/abs/1411.5672) treats propagation of linearized
  Regge data under sequences of simplicial moves and identifies lattice
  gravitons only after a gauge/curvature analysis.
- [Barrett et al.](https://arxiv.org/abs/gr-qc/9411008) and
  [De Felice--Fabri](https://arxiv.org/abs/gr-qc/0009093) implement multi-step
  Sorkin evolution for the dust-filled 600-cell, but do not publish the
  complete anisotropic three-slice Hessian below.

Thus the variational identity and its equivalence to canonical propagation
are **KNOWN**.  The proposed operator is not a new general formalism.

No located primary source prints the present `720`-position, fixed-mass,
staircase-carrier operator, its seven minimal blocks, or its calibrated twist
regularity.  External novelty remains **OPEN**; absence from a search is not
proof.

## KNOWN / CONTROL / OPEN

- **KNOWN:** a regular discrete Lagrangian gives equivalent symplectic
  one-step and three-point variational descriptions.
- **CONTROL:** both nonlinear slab solutions pass their complete internal
  equations and their thirty-component homogeneous canonical seam.
- **CONTROL:** the full maps `T_1`, `T_2` and `T_2 T_1` are already certified
  canonical on all `1440` phase dimensions in both schedules.
- **CONTROL:** all full pre-Legendre blocks are regular; this does not by
  itself prove the boundary-to-boundary twist blocks `B_i` are invertible.
- **OPEN:** invertibility and conditioning of every `B_i` in all seven
  sectors.
- **OPEN:** whether the Hessians reconstructed independently from the
  canonical blocks satisfy their adjoint identities inside propagated
  errors.
- **OPEN:** whether the three-slice recurrence reproduces the committed
  product map inside propagated errors.
- **OPEN:** schedule robustness of the Jacobi coefficients.
- **OPEN:** any decomposition into temporal and spatial parts.
- **OPEN:** continuum tensor modes, dispersion, a limiting speed,
  refinement and nonlinear anisotropic evolution.

## Framing attack

Constructing a three-slice recurrence is not yet deriving the wave equation.
Any regular pair of discrete canonical maps permits a position recurrence;
its existence alone has no gravitational content.  The scientific content
at this gate is narrower:

1. whether the action-selected boundary twist exists without a hidden
   singular direction;
2. whether the recurrence is genuinely the Hessian of the same two-slab
   action and not an algebraic manipulation inconsistent with the seam;
3. whether it survives the two schedule constructions.

The normalized `P,Q` depend on position coordinates, while the unnormalized
`K_-,K_0,K_+` retain the variational Hessian structure.  Neither spectrum is
a frequency spectrum on an evolving background.

A comparison with an intrinsic spatial operator must occur only after this
target-free operator is committed.  Even a later spectral correlation would
be a **PATTERN** unless a common symmetry decomposition selects the pairing.
A physical `c` additionally needs a derived clock normalization and
refinement; a ratio of arbitrary matrix norms is forbidden.

