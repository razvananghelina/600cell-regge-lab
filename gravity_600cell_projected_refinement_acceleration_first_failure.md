# Preserved first blind failure: the lapse leading equation is quadratic

Date: 2026-08-17

Protocol commit: `23157e2`  
Unit-normalization correction: `05d5685`  
Registered first implementation: `f474463`

First artifact:

```text
reproducible/gravity_600cell_projected_refinement_acceleration_blind.json
SHA-256 3b298170728c57fcdba26abcd31b74aa4dde13082ece1ba3d4be29047e3d3218
```

## Outcome

The first targeted run exited `7/9` with

```text
PROJECTED_REGGE_ACCELERATION_COEFFICIENTS_OPEN
```

No refined coefficient was compared with the continuum target.  The action,
carrier, volume-radius, exact coarse-action, Lorentzian branch, relabelling,
seam extrapolation and exact coarse coefficient controls passed.  The common
failure was the preregistered claim that the leading lapse residual is affine
in the acceleration coefficient.

## Why the failed hypothesis is mathematically false

The seam residual has an affine leading limit and produced the stable coarse
coefficient

```text
-0.539489745315062
```

against the exact calibrated value

```text
-0.5394897340206755...
```

But the coarse extrapolated lapse residuals were

```text
f(0)  = -0.000117645589...
f(-1) = 12.898858398...
f(-2) = 81.784388764...
```

They are plainly not affine.  This is not a physical disagreement between
the constraint and seam equations.  At fixed endpoint difference, the lapse
derivative sees the term

```text
(s_plus-s_minus)^2/rho.
```

With `s_plus-s_minus=O(a*eta^2)` and `rho=eta^2`, the leading lapse equation
is quadratic in `a`.  The static mass identity supplies one exact root
`a=0`; the nonzero root is the dynamic branch.  The first protocol had
incorrectly tried to infer an affine root from the static and one nonstatic
value, which necessarily returned a number near zero and failed its held-out
check.

## Numerical derivative issue exposed at the same time

The complex analytic-branch derivative is accurate for the seam momenta, but
is poorly conditioned for the lapse residual after the exact `O(eta)`
curvature/dust cancellation.  On the known level-zero control, a symmetric
real-log derivative with Richardson extrapolation at steps `2e-3` and `1e-3`
recovers the nonzero quadratic root stably, whereas the complex step does
not.  This was established using only the exact coarse calibration, before
using any refined target comparison.

## Frozen correction required before a rerun

The correction must:

1. retain the affine seam estimator unchanged;
2. add the held-out sentinel `a=-3`;
3. derive the lapse quadratic from `a=-1,-2`, retain the static root `a=0`,
   and compare its nonzero root with the seam root;
4. test the predicted `a=-3` lapse value;
5. use the preregistered real-log derivative and repeat step for the lapse;
6. rerun all four carriers without inspecting continuum distances.

This is a **DERIVED METHODOLOGICAL NEGATIVE** against the first lapse
estimator, not against refined Regge dynamics.  The failed artifact remains
in git history and cannot be silently replaced.
