# Preregistration: exact collective lapse path

Date: 2026-08-13

Upstream Schur result commit: `dc927a5`

Prior-art update: `0882934`

Status: **frozen before evaluating the arbitrary-precision action on the
exact collective path**.

The existing `FIVE_STIFF` label is not changed retroactively.  This protocol
tests a different, exactly specified direction because the preceding result
proved that its double-precision Schur lift omitted a relevant lift-error
term.

## 1. Frozen path

Keep the unrounded published `l0`, dust mass, old boundary and final boundary
fixed.  For a real path coordinate `t`, set

```text
rho(t) = tau^2 exp(t),
q(t)   = l0^2-rho(t),
```

on all five pole orbits and all thirty staircase-diagonal orbits,
respectively.  This is the published time-symmetric relation itself, not a
direction obtained from a numerical Hessian.

The base tangent in logarithmic internal coordinates is

```text
w_exact = (-tau^2/(l0^2-tau^2) repeated 30 times,
            1                         repeated  5 times).
```

No component is fitted to a residual or spectrum.

## 2. Frozen action and residual evaluations

Use the independently certified complete action-only evaluator at 100 decimal
digits.  Evaluate the action at exactly

```text
t in {-0.1, -0.03, -0.01, -0.003, -0.001,
       0,
       0.001, 0.003, 0.01, 0.03, 0.1}.
```

Report `S(t)-S(0)`, all imaginary parts and the maximum normalized action
difference.  The frozen constancy threshold is `1e-50` with denominator
`max(1,abs(S(t)),abs(S(0)))`.

At the same eleven points, evaluate the certified analytic complete-action
gradient and report all 35 local per-edge residuals.  Require every
representative simplex to remain Lorentzian and off branch boundaries.  The
frozen stationarity threshold is maximum absolute residual `1e-7`.

## 3. Frozen direct tangent curvature

Using the exact path, compute first and second centered differences at

```text
h = 1e-2, 5e-3, 2.5e-3
```

and the same two-stage Richardson/extrapolation rules as in the Schur
correction.  Report the first derivative, second derivative and step
convergence.  The frozen exact-path null threshold for the extrapolated
second derivative is `1e-40` in the action normalization.

This threshold tests computational constancy at 100 decimals; it is not an
interval proof of an exact identity.

## 4. Frozen four-mode control

Do not recalculate or select the four relative modes.  Load their fifteen
direction, 80-decimal reconstruction from commit `dc927a5` and require:

```text
min(abs(relative eigenvalues)) > 100 * epsilon_5
```

for both parities.  Record the values unchanged.

## 5. Frozen outcome labels

Assign exactly one label per parity.

### ONE_COLLECTIVE_LAPSE_NULL_FOUR_PSEUDOCONSTRAINT_STIFF

- all eleven full residual evaluations pass `1e-7`;
- all eleven path geometries pass the branch gates;
- the 100-decimal action is constant to `1e-50`;
- extrapolated exact-path curvature is at most `1e-40`;
- all four recorded relative modes remain separated by `100 epsilon_5`.

### TANGENTIAL_NULL_NOT_SOLUTION_FAMILY

The action/path-curvature gates pass but at least one transverse full residual
fails.

### COLLECTIVE_DIRECTION_NOT_NULL

Either the action constancy or exact-path curvature gate fails by more than
the frozen threshold while numerical evaluation remains otherwise resolved.

### NUMERICALLY_UNRESOLVED

Branch failure, precision contamination or any unclassified case.

The physical outcome is not a verifier PASS target.

## 6. Claim boundary

A positive `1+4` outcome would establish, computationally on this restricted
carrier, that:

- the chosen `tau` is gauge/lapse data rather than a selected physical tick;
- fixing five lapse variables leaves a regular thirty-equation evolution
  block;
- four relative phase-lapse combinations are small but resolved
  curvature-induced pseudo-constraints;
- one collective pole equation is redundant along the time-symmetric family.

It would **not** make the four relative modes particles, derive physical time,
prove a continuum constraint algebra, cover the 840-edge carrier, or produce
multi-tick dynamics.  Lapse nullity and pseudo-constraints are known Regge
phenomena; external novelty of this explicit 600-cell realization remains
**OPEN**.
