# Protocol: direct adversarial audit of shifted pseudo-longitudinal persistence

Date: 2026-08-19

The target-disclosed primary result is frozen as
`SHIFTED_PSEUDOLONGITUDINAL_DEFECT_PERSISTS`.  This protocol is committed
before any direct shifted residual is reconstructed.

## Purpose

Attack the near equality of current and shifted dimensionless residuals.  The
primary calculation used the committed shifted centered binary64 archive and
column-pivoted QR.  The adversarial calculation must instead rebuild both
shifted slabs from local 4-simplex action Hessians and use exact golden-ratio
geometry plus SVD.  Agreement is meaningful only if no centered `M,V` archive
or primary residual is used to construct the direct census.

## Frozen provenance

| input | SHA-256 |
|---|---|
| primary shifted verifier | `e4c5bcc18007c1c0ba7fbd38e29dffcc33a526fd790dbfcba8defe2ae44b7ab2` |
| primary shifted artifact | `0480f5d49d24e0f5d8e4e95f0cf62b7d0d9242459ed2b8f6d8e835ecd6e103a7` |
| accepted shifted direct-rank verifier | `1b54cd25899037fc66c2b58e01ef3bac267c6ebf2c6917d2a05ac4ac0feed1c5` |
| accepted shifted direct-rank artifact | `86b53f228d6cfa7326a677d881463f1b849e76bc6c9ac2b0e8aa6fd427042944` |
| prior-art gate | `740eefaee14ea3ff634f8cff237041cecd675c4ceaf3d5be6ccdb9a3778a57ef` |
| primary protocol | `f09c7450f238774b48425674e5b0800d6ef8ea2fee9f45b82794a0edef0a2375` |
| current direct action-York verifier | `73d852d58b21a9a15306a565d5cf4fb998b159fadb82830739ab0996ac07270e` |
| first accepted tick | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |
| second accepted tick | `936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70` |
| third accepted tick | `ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0` |
| binary-orbit verifier | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| anisotropic-rank source | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |
| boundary-tangent source | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| conformal source | `d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4` |
| 600-cell implementation | `ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f` |

## Independent construction

For both parities:

1. Read only the accepted first, second and third homogeneous states needed
   to define slabs 2 and 3.
2. Reconstruct every local Hessian pattern at 100 decimal digits for steps
   `1e-20, 1e-15, 3e-20, 3e-15`, retaining the accepted action-Hessian
   Hermitian projection and schedule-envelope control.
3. Assemble the global representative kernels and solve the two direct
   Legendre maps.  Form direct shifted

   ```text
   M = (K_minus + K_plus)/2,
   V = K_minus + K_zero + K_plus.
   ```

4. Construct the 120 exact normalized golden-ratio vertices independently of
   `commons.cell600`, then require exact agreement with its frozen labels and
   all 720 edges.
5. Construct the sector-compressed incidence `C` and tangential
   vertex-displacement matrix `D` directly from these vertices.  Use sectors
   4 and 5 only.
6. Use full SVD—not the primary QR path—to derive the conformal complement
   `W` and the 15-dimensional longitudinal image `L`.
7. Compute the same dimensionless residuals

   ```text
   rho_span = ||(1-P_BL) A L|| / ||A L||,
   rho_comm = ||(1-P_L) B^-1 A L|| / ||B^-1 A L||,
   ```

   and the augmented singular spectrum of `[B L,A L]`.

No shifted centered NPZ and no primary pseudo-longitudinal JSON may be opened
before all 16 direct cells have been calculated and classified.  The primary
JSON may then be opened only for a post-census comparison.

## Direct error and classification

For each cell, retain the direct `acb` midpoint and component re-enclosure
radii.  The local dimensionless floor is

```text
1000 eps 30 kappa
+ kappa * max(epsilon_M / max(1,||M||),
              epsilon_V / max(1,||V||)),
```

where `epsilon_M,V` are the existing conservative matrix-error bounds and
`kappa` is the maximum of the carrier, kinetic, `BL` and normalization
condition factors.  For each parity/sector family, the final error is the
maximum four-schedule variation from `operational_primary` plus the maximum
local floor.

The frozen labels are

```text
value <= 10 error     ZERO_CONSISTENT
value > 100 error     NONZERO_RESOLVED
otherwise             OPEN.
```

The augmented-rank threshold is `100 * max(final errors) * sigma_max`; only
rank equal to 15 or strictly greater than 15 is interpreted.  Both exact norm
inequalities relating span and inverse residuals must pass their numerical
floor.  All kinetic forms must be positive and all generalized stiffness
inertias must be `15 negative + 10 positive`.

## Frozen outcome hierarchy

1. `SHIFTED_PSEUDOLONGITUDINAL_DIRECT_CONTROL_FAILED` if provenance,
   state, branch, Hessian, Legendre, exact-geometry, carrier, denominator,
   conditioning or completeness controls fail.
2. `SHIFTED_PSEUDOLONGITUDINAL_DIRECT_PERSISTENCE_CONFIRMED` if both
   residuals are `NONZERO_RESOLVED` and augmented rank is greater than 15 in
   all 16 direct cells.
3. `SHIFTED_PSEUDOLONGITUDINAL_DIRECT_PERSISTENCE_REFUTED` if both residuals
   are `ZERO_CONSISTENT` and augmented rank is 15 in all 16 cells.
4. `SHIFTED_PSEUDOLONGITUDINAL_DIRECT_OPEN` otherwise.

## Interpretation firewall

Confirmation establishes only a mechanically independent two-tick finite
background result.  It does not make the near numerical equality a conserved
quantity, because only two times are known and the residual is invariant
under broad scalar rescalings.  Curvature scaling, refinement, continuum
recovery, physical instability, propagation and speed remain **OPEN**.

Only this targeted verifier and static registry guards may be run.  The full
suite remains excluded by the user's standing instruction.

