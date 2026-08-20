# First complete-carrier/action intersection result: numerically open

Date: 2026-08-20

## Outcome

The corrected targeted verifier passed `13/13` with

```text
FULL_SCALE_STRUT_CANONICAL_NUMERICALLY_OPEN.
```

**OPEN, not a physical result.**  Every provenance, geometry, pole-identity,
input-rank, reconstruction and corruption control passed, but the frozen
binary64 calibration cannot resolve any of the 14 parity/sector
intersections.

## Provenance

- prior-art gate: `c2d0e83`;
- target-blind protocol: `b621736`;
- registered verifier: `8e20920`;
- frozen first classifier failure: `61fc3cc`;
- classifier correction protocol: `cf0a296`;
- classifier-only repair: `9522e12`;
- frozen OPEN artifact: `91b337b`;
- artifact SHA-256:
  `b29cc33a9effeb2087fb6133359ee747d100d203778586372a7ceeebc2e4f070`.

No full suite was run.

## What the calculation says

The independent analytic and stored carrier reconstructions agree with
maximum coefficient discrepancy `2.44e-65` and binary matrix discrepancy
zero.  The action-response source reproduces byte-identically.  In every
minimal sector the literal pole blocks and the ranks of the input graphs are
correct.

For the six non-homogeneous sectors, the smallest operational scaled `D`
singular values lie around

```text
0.96e-7 ... 1.27e-7,
```

or only `11.88 ... 15.65` times the frozen `epsilon ~= 0.81e-8`.  They are
above the `10 epsilon` zero threshold but below the `100 epsilon` nonzero
threshold.  Exactly the `5d` weaker singular directions in each such sector
remain open.

The homogeneous trivial sector has one reduced singular value around
`3.86e-16`, while four further values remain in the same `1e-7` open band.
The reduced matrix operationally reports one zero but the more redundant
joined matrix reports two.  Exact pole algebra says their true nullities must
agree, so this disagreement is a conditioning warning, not evidence for two
modes.

Both parities reproduce the same qualitative pattern.  That agreement is a
control, not a rank proof.

## Why the calibration is insufficient

About `8.10e-9` of each `epsilon` comes from applying the global carrier's
`5.73e-9` high-precision **spectral discrepancy** as if it were a matrix
perturbation bound.  Flint lift radii and action derivative-step differences
are only around `1e-84` and `1e-30` after scaling.  The proxy was frozen in
the protocol and therefore correctly forces an OPEN verdict, but it is too
coarse to decide the new intersection.

## Post-result literature check

- Dittrich, Freidel and Speziale, [*Linearized dynamics from the 4-simplex
  Regge action*](https://arxiv.org/abs/0707.4513), explicitly connect Regge
  Hessian zero modes with remnant diffeomorphism symmetry and continuum
  graviton analysis.
- Hoehn, [*Canonical linearized Regge Calculus*](https://arxiv.org/abs/1411.5672),
  requires gauge analysis before calling curvature modes lattice gravitons.
- Dittrich, Kaminski and Steinhaus,
  [arXiv:1404.5288](https://arxiv.org/abs/1404.5288), analyze when Regge
  Hessians become singular and emphasize the geometric/background
  dependence of such claims.

These sources support the need to distinguish an exact Hessian null mode
from numerical near-nullity.  They do not resolve this 600-cell carrier.
External novelty remains **OPEN**.

## Next falsifiable step

Rebuild the projected carrier directly in multiprecision and retain the
canonical response as Flint balls rather than binary midpoints.  Use
coordinate-scaled Gram determinants/minors to certify full rank where
possible.  The target is disclosed: non-homogeneous directions look nonzero
and one homogeneous direction looks null.  Neither is accepted unless the
high-precision audit can falsify the opposite.

