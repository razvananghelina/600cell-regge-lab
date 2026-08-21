# Primary result: classical tick scale covariance

Date: 2026-08-21

## Provenance

- prior-art gate commit: `3275575`;
- frozen protocol commit: `fd1f8a8`;
- registered implementation before first evaluation: `127a4af`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_tick_scale_covariance.py`;
- artifact:
  `reproducible/gravity_600cell_tick_scale_covariance.json`;
- artifact SHA-256:
  `8eeb9e76dd23672a8f055ce4119ab60d1cdbea8932d03a8cb021f945f8816f0e`.

Only the targeted verifier was run.  It returned **12/12**.  The full suite
and all deferred nonlinear censuses were not run.

## Mechanical outcome

```text
TICK_SCALE_COVARIANCE_PRIMARY_CONFIRMED
```

For both independently derived order-24 staircase parities, at the frozen
off-shell nonhomogeneous state and at `alpha=3/5` and `alpha=7/4`, the
100-decimal complete action obeys

```text
S(alpha^2 q, alpha M) = alpha^2 S(q,M),
g_log(alpha^2 q, alpha M) = alpha^2 g_log(q,M)
```

for all 95 logarithmic squared-length derivatives.  The maximum action errors
range from `4.48e-98` to `2.12e-93`; the maximum component errors range from
`5.83e-99` to `1.70e-95`, against the frozen `1e-65` gate.

The off-shell support gate is nonvacuous: all 35 internal and all 60 boundary
components exceed `1e-20` for both parities.

## Exact algebra and mass rules

SymPy certifies the scale degrees of the triangle area polynomial, generic
4-simplex Gram determinant, dihedral sine and cosine, Regge hinge term and
dust term.  It also certifies length degree one for both masses currently used
by the project:

```text
M_coarse = (90/pi) epsilon3 L,
mu_r     = K_r/(8*pi).
```

The second certificate is tied to the accepted refined curvature-mass artifact,
not inferred from a numerical coincidence.

## Falsification control

Holding the mass fixed while scaling every squared length destroys the
identity, as it must.  The normalized action defects are `0.62768` and `1.0`;
the pole-gradient defects range from `0.0026153` to `0.014303`.  Thus the
positive result is sensitive to the matter scaling hypothesis and cannot be
explained by a test that always returns covariance.

## Provisional interpretation

- **DERIVED EXACT / PRIMARY:** the stated action is homogeneous of length
  degree two when its geometrized dust masses scale with the geometry.
- **DERIVED COMPUTATIONAL / PRIMARY:** the complete orbit implementation
  realizes that identity off shell on both parities.
- **STRUCTURAL:** zeros of its internal equations and compatible canonical
  data therefore transform in continuous global scale families.
- **NOT YET ACCEPTED:** the absolute-tick no-go awaits the separately frozen
  direct-2400-simplex adversarial audit.
- **NOT EXCLUDED:** selection of `tau/L`, `tau_next/tau0`, relational dust
  time, or a duration relative to an externally fixed dimensional anchor.
- **NOT DERIVED:** seconds, `c`, `G`, `hbar`, Planck time or a quantum scale.

The result is conditional.  A fixed mass, cosmological constant, boundary
length or quantum scale breaks the stated hypotheses by supplying a scale.

