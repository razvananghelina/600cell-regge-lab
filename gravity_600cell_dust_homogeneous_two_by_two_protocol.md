# Preregistered protocol: homogeneous `2 x 2` tangent and curvature kernel

Date: 2026-08-17

Prior-art gate commit: `24eed99`.

This is a target-disclosed confirmatory test.  The previously observed raw
Rayleigh multiplier near `-1` is disclosed in advance.  No result from the
new `2 x 2` reduction has been inspected before this protocol commit.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_homogeneous_curvature_kernel.json` | `b55887ff3905afd94e86821852d58f0d60c227b52dfbd945044874bfe87540e9` |
| `verify_gravity_600cell_dust_homogeneous_curvature_kernel.py` | `43837b4d97fcf21cc6de9e4debea0c22bc827d5186d0c75ca07dfe5c799e1a15` |
| `gravity_600cell_dust_internal_curvature_response.json` | `95b6edd8e21ad20a0db97a7c8e7027db7da6547b2b994ad1eb595cf2307f29dc` |
| `verify_gravity_600cell_dust_internal_curvature_response.py` | `276982879fae5f8fa735f27a6fa30bfe965dc3e41c169d8a229a61c23511ae66` |
| `gravity_600cell_dust_full_boundary_tangent.npz` | `816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `verify_gravity_600cell_dust_hyperbolic_lapse_alignment.py` | `e461296a965c9b80fb89fae5660ce642858f3d3dfa0b24ccdecc2aced53c7047` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| `gravity_600cell_dust_homothetic_canonical_lapse.json` | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |

Both schedule parities and the four already frozen derivative variants are
mandatory.  Arithmetic uses 100 decimal digits and 80-decimal Flint balls.
Only the trivial sector is reconstructed.

## Fixed carrier

Let

```text
u_q = (1_30,0_30)/sqrt(30),
u_p = (0_30,1_30)/sqrt(30),
U   = [u_q u_p],
P   = U U*.
```

This basis is fixed before computation.  It is not rotated toward the
observed kernel.  For each parity and derivative variant reconstruct the
same ball tangent `T : C^60 -> C^60` and curvature response
`F : C^60 -> C^160` used by the two passing parent verifiers.

Form

```text
A = U* T U                         (2 x 2 compression),
L = (I-P) T U                      (58 x 2 plane leakage),
B = F U                            (160 x 2 curvature response).
```

All midpoints and propagated ball radii are retained before any binary64
diagnostic.  The `2 x 2` eigensystem and the `B*B` kernel calculation use
high-precision arithmetic, not the ill-conditioned 60-dimensional
eigensystem.

## Frozen calibration

For any reported norm or principal-line distance, its uncertainty is the sum
of:

1. operational-primary versus operational-shadow change;
2. validation-primary versus validation-shadow change;
3. operational-primary versus validation-primary change;
4. the maximum propagated Flint radius over all four variants;
5. a dimension-scaled high-precision/binary serialization floor;
6. `1e-70`.

For a nonnegative diagnostic `x` with calibrated uncertainty `epsilon`, use

```text
ZERO / IDENTIFIED       x <= 10 epsilon,
NONZERO / SEPARATED     x > 100 epsilon,
NUMERICALLY_OPEN        otherwise.
```

No raw small decimal overrides this rule.  Uncertainty greater than or equal
to `1e-2` for a bounded angular distance forces `NUMERICALLY_OPEN`.

## Controls and tests

For each parity:

1. reproduce the exact frozen provenance, the 43 imported geometry controls,
   the unique trivial sector, the 65-dimensional pre-Legendre determinant
   excluding zero, and parent rank `59/nullity 1`;
2. require `U*U=I_2` at the arithmetic floor;
3. classify `||L||_2`.  `U` is called invariant only when the leakage is
   `ZERO`; it is called non-invariant only when leakage is `NONZERO`;
4. classify `B` as rank one only when its larger singular value is resolved
   nonzero and its smaller singular value is zero under the frozen
   calibration;
5. obtain the unit kernel vector `k in C^2` from the smaller eigenline of
   `B*B`, with angular uncertainty divided by its resolved spectral gap;
6. test the *full* eigenline residual

   ```text
   r_K = ||(I - U k k* U*) T U k||_2,
   mu  = k* U* T U k.
   ```

   This test remains meaningful even if the rest of `U` leaks;
7. diagonalize `A` independently using its quadratic characteristic
   polynomial and report both eigenvalues, determinant, trace, symplectic
   defect `A* J A - J`, and their calibrated uncertainties;
8. classify `|mu+1|` against its propagated uncertainty.  The labels are
   `EXACT_MINUS_ONE_WITHIN_ERROR`, `RESOLVED_NOT_MINUS_ONE`, or
   `NUMERICALLY_OPEN`.  Nearness alone is not equality.

Across schedules, compare the two kernel vectors in their common `(q,p)`
basis and compare the two `A` matrices.  The earlier literal orbit-set
permutation must be reproduced independently; because it is a permutation,
it fixes both normalized uniform vectors exactly.  A resolved schedule
difference outranks every positive interpretation.

## Frozen outcome hierarchy

1. `HOMOGENEOUS_2X2_CONTROL_FAILED` for any provenance, reconstruction,
   determinant, rank or calibration-control failure.
2. `HOMOGENEOUS_2X2_SCHEDULE_DEPENDENT` if either kernel line or compressed
   matrix is resolved different between schedules.
3. `HOMOGENEOUS_2X2_CURVATURE_KERNEL_OPEN` if `B` is not resolved rank one in
   both schedules.
4. `HOMOGENEOUS_2X2_KERNEL_NOT_EIGENLINE` if the full residual `r_K` is
   resolved nonzero in either schedule.
5. `HOMOGENEOUS_2X2_EIGENLINE_PLANE_LEAKS` if the kernel is an identified
   eigenline but the full homogeneous plane has resolved leakage.
6. `HOMOGENEOUS_2X2_EIGENLINE_PLANE_OPEN` if the kernel is an identified
   eigenline but invariance of the full plane is numerically open.
7. `HOMOGENEOUS_2X2_INVARIANT_EIGENLINE` only if both the kernel line and the
   whole plane are invariant in both schedules.
8. `HOMOGENEOUS_2X2_NUMERICALLY_OPEN` otherwise.

The separate `mu=-1` label does not choose the main outcome.

## Explicit exclusions

- no 60-dimensional eigenvector is used to define or stabilize `k`;
- no fitted basis rotation, threshold or schedule conjugacy;
- no inference that an eigenline is gauge, lapse, time or a graviton;
- no nonlinear continuation, proper-time normalization, dispersion or
  continuum comparison;
- no full-suite run.
