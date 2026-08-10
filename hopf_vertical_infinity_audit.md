# Hopf vertical infinity audit

Date: 2026-08-10

## Question and complete hypotheses

Let the unit round three-sphere carry its standard Hopf fibration

`S1 -> S3 -> S2`.

On complex scalar functions, write the non-negative Laplacian as

`Delta = Delta_H + Delta_V`,

where `Delta_V=-K_z^2` differentiates along the circle fibers and
`Delta_H=-(K_x^2+K_y^2)` is the horizontal sub-Laplacian.  This is a
continuum statement about the round Hopf geometry.  It is not yet a statement
about Lorentzian time, particles, or the finite 600-cell.

Robert O. Bauer derives the same simultaneous spectrum in *Analysis of the
horizontal Laplacian for the Hopf fibration*, Forum Mathematicum 17 (2005),
903--920, arXiv:math/0209127.

## Exact simultaneous spectrum

The scalar harmonics are indexed by

`N=1,2,...`,  `k=-(N-1),-(N-3),...,N-1`.

For every fixed `(N,k)`, the weight space has complex multiplicity `N` and

`lambda_full = N^2-1`,

`lambda_V = k^2`,

`lambda_H = N^2-1-k^2`.

There are `N` allowed weights and hence total multiplicity `N^2`, as required
for the degree-`N-1` scalar harmonics of `S3`.

## Does the infinity exist?

Yes.  The vertical eigenvalue is zero whenever `k=0`, which occurs for every
odd `N`.  The cumulative kernel multiplicity through `N<=L` is

`sum_(odd N<=L) N = ceil(L/2)^2`.

It diverges.  Geometrically these are precisely pullbacks of functions on the
Hopf base: every `f on S2` gives `f compose pi` constant along each fiber.
Thus `ker Delta_V` contains `L2(S2)` and is infinite-dimensional.

The problem is stronger than an infinite kernel.  The positive vertical
eigenvalue `1` occurs at `k=+1` and `k=-1` for every even `N`, so it too has
infinite multiplicity.  Therefore `Delta_V` alone has no compact resolvent and
is not an elliptic spectral geometry on `S3`.

This infinity is **DERIVED**, not a numerical extrapolation.

## Is the infinity physically needed?

No.  It is not a prediction of infinitely many particles.  It says that the
vertical derivative alone ignores all variation on the two-dimensional base.
It is an incomplete kinetic operator.

The horizontal operator has only the constant zero mode.  Indeed, for an
allowed weight,

`lambda_H >= N^2-1-(N-1)^2 = 2(N-1)`,

so it vanishes only for `(N,k)=(1,0)`.  Its first positive eigenvalue is `2`,
on the four-dimensional `N=2`, `k=+-1` sector.

For any fixed `r>0`, the combined anisotropic operator

`Delta_r = Delta_H + r Delta_V`

has eigenvalues

`N^2-1-k^2+r k^2`.

Both summands are non-negative and their kernels intersect only in the
constant function.  Hence `ker Delta_r=C`, and its eigenvalues escape to
infinity.  The round metric selects `r=1`, giving the ordinary scalar
Laplacian.  Other positive `r` values are Berger-type anisotropies and require
an additional geometric selector.

Therefore:

- **DERIVED:** the separated vertical operator has an infinite kernel and
  infinite positive multiplicities;
- **DERIVED:** every positive combined operator removes that degeneracy down
  to the constant mode;
- **STRUCTURAL:** treating the two separated first gaps as two physical
  propagation speeds was the wrong observable;
- **OPEN:** no internal principle currently selects a non-round coefficient
  `r`, and inserting `r=a_1=5` would be circular.

The formal continuum component-gap quotient is `2/1=2`, not `5`, but it is
not promoted to a physical observable: the denominator belongs to a
non-elliptic operator with infinite multiplicity.

## What the current finite refinement does

The preregistered 600-cell split has exact finite kernel counts

```text
                           coarse     first barycentric level
fiber kernel                  12                 12
cross kernel                   1                  1
```

The twelve coarse fiber zero modes are one constant per discrete decagon.
Their failure to grow at the first barycentric level has a simple construction
reason: that refinement freezes one parentwise constant tensor and does not
refine the Hopf quotient or its twelve fibers.  It is a valid Galerkin
refinement of that piecewise-tensor PDE, but it is not yet a certified
fiber-adapted approximation to the smooth Hopf vector field.

This does not undo the earlier negative: the preregistered local-tensor
first-gap ratio really changes from `5` to `5.3388401713`.  It changes the
next task.  A continuum-faithful test must define the smooth Hopf field at new
geometric nodes and preregister a low-energy counting test.  Exact kernel
dimension alone is too brittle for non-fiber-aligned finite elements; the
relevant signal is a growing near-zero vertical band together with convergence
of the positive combined operator.

## Decision

**DERIVED CORRECTION:** the infinity exists, but only in the deliberately
separated vertical operator.  We do not need or want it in the physical
kinetic operator.  Its role is diagnostic: a proposed Hopf refinement must
explain how base modes are represented, while the combined operator must keep
only the constant zero mode.

The next construction is therefore fixed conceptually before computation:

1. use the smooth Hopf vector field, rather than a frozen parent-edge tensor;
2. refine/project the geometry on `S3`;
3. measure the vertical near-zero counting function and the combined spectrum;
4. do not compare any new scalar with `a_1=5` unless an independent selector
   first singles it out.
