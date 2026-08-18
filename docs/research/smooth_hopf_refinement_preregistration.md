# Preregistration: smooth orthogonal Hopf refinement

Date: 2026-08-10

## Why a new construction is required

The committed edge-supported split is exact combinatorics, but it is not the
continuum Hopf vertical/horizontal decomposition.  On each parent tetrahedron
its tensors have tangent eigenvalues

`spec(Q_f)=(1/2,0,0)`, `spec(Q_c)=(1/2,1,1)`.

In particular `Q_c` is rank three.  It represents the other five edge
contributions and still penalizes the fiber direction.  The true orthogonal
Hopf tensors have ranks one and two.

This preregistration fixes a continuum-faithful alternative before its
refined mode data are computed.  It contains no bootstrap integer, speed, or
target ratio.

## Geometric definition fixed before execution

Choose the same discrete Hopf fibration used by the blind edge census.  If its
right order-ten generator is

`g=cos(pi/5)+u sin(pi/5)`, `|u|=1`,

then the corresponding smooth unit Hopf vector field on the round sphere is

`X(q)=q*u`.

The invariant Hopf map is

`pi(q)=q*u*q^-1 in S2`.

At every mesh level:

1. barycentric vertices are radially projected back to the unit `S3`;
2. every child is the chordal tetrahedron joining its four projected nodes;
3. at the normalized spherical centroid of a tetrahedron, project `X` into
   the affine tangent three-plane and normalize it to `X_t`;
4. set

   `P_V=X_t tensor X_t`, `P_H=P_tangent-P_V`;

5. assemble consistent `P1` mass and stiffness matrices with one centroid
   quadrature point per tetrahedron.

The quadrature is fixed, not chosen after spectra.  Since basis gradients are
constant and the tensor field is smooth, this is a convergent low-order rule.
It also enforces the exact finite identity

`K_V+K_H=K_full`

on each mesh.

One representative fibration is sufficient for this test: the six discrete
fibrations are related by 600-cell symmetry.  No claim that geometry selects
one physical fibration is made.

## Continuum calibration modes fixed before execution

Two canonical function spaces avoid an arbitrary near-zero cutoff.

### Fiber-charged first harmonics

The four ambient coordinate functions of `q in S3` span the `(N=2,k=+-1)`
sector.  Their exact continuum eigenvalues are

`(lambda_V,lambda_H,lambda_full)=(1,2,3)`.

### Fiber-invariant base harmonics

The three components of `pi(q)` are pullbacks of the first coordinate
harmonics on the Hopf base.  They span `(N=3,k=0)` and have exact eigenvalues

`(lambda_V,lambda_H,lambda_full)=(0,8,8)`.

For each space the verifier records the generalized Ritz eigenvalues of the
three assembled forms on both the coarse and first projected-barycentric
meshes.  It also records low unconstrained spectra, but no target comparison
is made in the blind artifact.

## Preregistered comparison criteria

The later comparison will use only the following gates.

1. **Algebraic gate:** local projectors have ranks `1/2`; their sum and the
   assembled split identity hold to numerical precision.
2. **Geometric gate:** the maximum chord length decreases under refinement and
   the projected Hopf field never degenerates in an element tangent plane.
3. **Charged-mode gate:** the maximum absolute errors from `(1,2,3)` decrease
   from coarse to fine.
4. **Base-mode gate:** the maximum vertical Ritz value decreases toward zero,
   while the horizontal/full errors from `8` decrease.
5. **Combined-operator gate:** the low full spectrum has only the constant
   zero mode.

Failure of any convergence-direction gate is a negative for this particular
projected `P1` implementation, not a theorem against the smooth Hopf
decomposition.  Passing all gates would validate the discretization, not
select a physical anisotropy or recover `a_1=5`.

## Explicit exclusions

- No separated gap ratio is an acceptance target.
- No coefficient `r` in `K_H+r K_V` is varied or fitted.
- No Lorentzian interpretation of the compact Hopf circle is made.
- No result from the earlier `5 -> 5.33884` comparison is used to choose a
  cutoff, quadrature, or mode space.

## Post-result scope correction (not part of the preregistered gates)

The phrase "convergent low-order rule" above was too broad.  The fixed
quadrature is consistent on the two stated meshes, but a convergence theorem
for an infinite tower additionally requires shape regularity.  Repeated
barycentric refinement is already known in this repository to make that
nontrivial.  The later result can therefore establish a successful first-level
calibration, not infinite-tower convergence.  No acceptance gate or numerical
definition was changed by this correction.
