# Blind dynamic dust tangent census

Date: 2026-08-17

## Status ledger

- **DERIVED:** the two accepted first-tick schedules define calibrated, full-rank
  canonical tangent maps on the order-24 invariant quotient after the 35 internal
  variables are eliminated.
- **DERIVED:** both 60-dimensional maps satisfy the symplectic test at the
  preregistered tolerance.
- **DERIVED:** the global scale plane and the 58-dimensional shape phase space
  are invariant at the preregistered tolerance for both schedules.
- **STRUCTURAL:** this is the order-24 invariant quotient with 30 boundary edge
  orbits.  It is not the full 720-edge carrier and it does not establish a
  graviton interpretation.
- **OPEN:** the frozen mechanical classifier says that the two shape spectra are
  schedule-dependent, but the separation was obtained with binary64
  eigendecomposition of nonnormal matrices.  A high-precision, basis-independent
  audit is required before the label is given physical meaning.
- **OPEN:** continuum matching, a limiting speed, refinement, long-time
  stability, and emergent time were not tested.

## Frozen provenance

The calculation was performed only after the following commits:

- prior-art gate: `25722d9`
- blind protocol: `0bceb9b`
- registered verifier implementation, before its first run: `b79b4a3`

The protocol explicitly forbade parsing a continuum spectrum or a speed target.
The artifact records `continuum_target_parsed = false`,
`speed_target_parsed = false`, and `full_720_edge_carrier = false`.

Verifier:
`reproducible/verify_gravity_600cell_dust_dynamic_tangent.py`

Artifact:
`reproducible/gravity_600cell_dust_dynamic_tangent.json`

Artifact SHA-256:
`1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5`

Only the targeted verifier was run.  It passed `12/12`; the full suite was not
run, following the explicit instruction for the current work.

## Blind results

For both time-order schedules, the 95-variable dynamic Hessian reproduced the
accepted nonlinear tick, passed the Lorentzian-domain checks for all derivative
evaluations, and passed the calibrated entry, operator-proxy, and validation-step
comparisons.  The reduced pre-Legendre Jacobian has resolved rank 65 and no
calibrated zero singular values.

The induced 60 by 60 canonical maps have:

| quantity | even schedule | odd schedule |
|---|---:|---:|
| symplectic defect | within `6.876e-18` tolerance | within `6.876e-18` tolerance |
| scale/shape mixing | `5.2755e-30` | `5.2755e-30` |
| mixing tolerance | `5.2755e-19` | `5.2755e-19` |
| full spectral radius | `45.37453401114913` | `45.37453401114913` |
| shape eigenvector condition estimate | `9.1513e4` | `5.6427e5` |

The exact preregistered scale/shape test therefore returns
`SCALE_SHAPE_INVARIANT` for both schedules.  This is evidence that the dynamic
linearization contains a separate 58-dimensional shape sector rather than only
the homogeneous Friedmann scale mode.

The frozen classifier returned
`DYNAMIC_SHAPE_TANGENT_SCHEDULE_DEPENDENT` because the optimally matched shape
eigenvalue distance is `6.8628065e-11`, above its stored eigenvalue uncertainty
`9.0749068e-14`.  The matched singular-value distance is `5.8207661e-11`, below
its stored uncertainty `8.9907519e-10`.

## Hostile numerical reading

The mechanical schedule verdict is not yet a scientific schedule-dependence
claim.  The two shape matrices are nonnormal: their reported eigenvector
condition estimates are approximately `9.15e4` and `5.64e5`.  Binary64 roundoff
amplified by those condition estimates can plausibly reach the scale of the
observed `6.86e-11` eigenvalue separation.  In contrast, the singular spectra do
not resolve a schedule difference under the frozen calibration.

This caveat does not alter the preregistered output after seeing it.  It creates
a separate correction mission: compare basis-independent spectral invariants
and, if needed, eigenvalues at high precision with an error rule frozen before
that calculation.  Until that audit is complete, schedule dependence is
**OPEN**, while full rank, symplecticity, and scale/shape invariance remain
**DERIVED** within the registered finite calculation.

## What this does not show

The radius `45.3745` is not by itself a physical instability: no physical norm,
constraint quotient beyond the stated reduction, continuum mode identification,
or refinement limit has been derived.  Likewise, discreteness, a spectral gap,
or degeneracies would not alone be new physics because compact continuum
three-spheres already have discrete spectra.  No claim about gravitational
waves, a physical tick, a causal cone, or the speed of light is made here.
