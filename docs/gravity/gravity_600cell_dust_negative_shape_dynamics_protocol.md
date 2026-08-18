# Target-disclosed protocol: complete negative-shape sector dynamics

Date: 2026-08-18

Status: **PREREGISTERED BEFORE ANY NEGATIVE-SUBSPACE LEAKAGE OR SECTOR
COMPANION WAS COMPUTED.**

The two sector labels and their `15+10` stiffness inertia are already known
from the preceding blind result.  This is therefore an explicitly
confirmatory/diagnostic protocol.  No companion eigenvalue, singular value,
invariance residual or desired stability count was inspected before it.

## Frozen provenance

```text
prior-art gate commit
  db1fa95

reproducible/gravity_600cell_dust_shape_stiffness.json
  03b1ad6bcc21af6481120ae00f04cbc06423f54ca5623cc5e0e2a251bd798868

reproducible/verify_gravity_600cell_dust_shape_stiffness.py
  d4f0a9a805910de37011ba70f407907daa2d11c650aeea22e571ab867282a44c

reproducible/gravity_600cell_dust_centered_jacobi.json
  fe0c2d231c2b7eaa8a96cc051de8b3a9b034e384589ab6411db81562af0d9b56

reproducible/gravity_600cell_dust_centered_jacobi.npz
  1077fb562abd4b16a9b5d664d5b7669e2ace0344022aa12bc071fcc4fd4691ef

reproducible/gravity_600cell_dust_full_anisotropic_legendre_rank.json
  7dc33fcebe8e2cb62be9bba5dfd1fca06fa176a06afe3717d2e9e866f67a7226

reproducible/verify_gravity_600cell_dust_conformal_supermetric.py
  d77dc8853826d9aecc4395fc4aae405d0505bbd644ec3a3229f640b2e980bcb4

reproducible/verify_gravity_600cell_dust_full_boundary_tangent.py
  c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571

reproducible/verify_gravity_global_regge_orbits.py
  ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf

commons/cell600.py
  ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f
```

## Complete hypotheses

1. the fixed regular 600-cell and its `720` logarithmic signed-squared spatial
   edge coordinates;
2. the first two accepted nonstationary fixed-total-mass dust-Regge slabs;
3. the literal adjacent-slice identification used by the centered Jacobi
   recurrence;
4. both staircase schedules and all four frozen derivative variants;
5. the canonical conformal incidence and action-relative shape complement;
6. minimal sectors `4` and `5`, selected only because the preceding committed
   result gives each the resolved inertia `15 negative + 10 positive`;
7. the committed midpoint/radius enclosures, `1000 eps n` arithmetic floors
   and `10/100` classification bands;
8. the complete Legendre rank `1560/1560`, so no exact fixed-carrier
   constraint quotient is applied;
9. no independent dust perturbation, physical time unit, refinement,
   continuum harmonic, polarization, speed or Planck target.

## Complete enumeration

The primary cells are

```text
2 schedules * 2 selected sectors * 4 variants = 16.
```

Each cell contains:

- two negative-subspace invariance tests (`Gamma_S`, `Omega_S`), for `32`
  total;
- one complete `50 x 50` sector companion, for `800` directly classified
  eigenvalue instances;
- if and only if both invariance tests pass, one canonical `30 x 30`
  negative-carrier companion, for `480` eigenvalue instances;
- complete ordered singular values for both available companions.

Both selected irreps have dimension one, so no additional irrep multiplicity
factor is required.

## Carrier reconstruction

Reconstruct `U`, `W`, the restricted Hermitian forms and all source errors
exactly as in the committed shape-stiffness verifier.  In each selected cell,
the inertia of

```text
A = -W*[(V+V*)/2]W
```

must reproduce `15 negative resolved + 10 positive resolved`, with no zero or
open entry.  Let `E` be the orthonormal basis of its negative eigenspace.

If the positive/negative gap is `g` and the complete restricted-form error is
`epsilon_A`, use

```text
eta_E = 2 epsilon_A / (g - 2 epsilon_A)
        + 1000 eps * 25
```

when `g>2 epsilon_A`; otherwise the negative carrier is open.

## Primary invariance tests

For `X` equal to `Gamma_S` and `Omega_S`, compute

```text
r_X = ||(I-EE*) X E||_2.
```

With restricted matrix error `epsilon_XS`, use

```text
epsilon_inv = epsilon_XS
              + 2 eta_E (||X||_2 + epsilon_XS)
              + 1000 eps * 25 * max(1,||X||_2).
```

Classify each residual:

```text
r <= 10 epsilon_inv       INVARIANT_CONSISTENT
r > 100 epsilon_inv       MIXING_RESOLVED
otherwise                 INVARIANCE_OPEN.
```

No projected negative companion is constructed for a cell unless both
operators are `INVARIANT_CONSISTENT`.

## Complete sector companion

For each cell define

```text
G = Gamma_S,  O = Omega_S,  R = I+G.
```

Classify `R` by its smallest singular value and complete error `epsilon_GS`:

```text
sigma_min(R) > 100 epsilon_GS   REGULAR_RESOLVED
sigma_min(R) <= 10 epsilon_GS   SINGULAR_CONSISTENT
otherwise                       REGULARITY_OPEN.
```

Only a resolved-regular cell receives

```text
D = R^-1,
L = -D(I-G),
H =  D(2I-O),
C = [0 I; L H].
```

Use the standard inverse perturbation envelope

```text
epsilon_D = ||D||^2 epsilon_GS / (1-||D|| epsilon_GS)
```

when the denominator is positive.  Propagate

```text
epsilon_L = epsilon_D (||I-G||+epsilon_GS)
            + ||D|| epsilon_GS + roundoff,

epsilon_H = epsilon_D (||2I-O||+epsilon_OS)
            + ||D|| epsilon_OS + roundoff,

epsilon_C = epsilon_L + epsilon_H + block roundoff.
```

All norms are spectral norms; roundoff is the frozen
`1000 eps n max(1,relevant norms)` envelope.

## Negative-carrier companion

Only after both invariance residuals pass, compress

```text
G_- = E* G E,
O_- = E* O E
```

and apply the same formula to obtain `C_-` on `E direct-sum E`.  Its error
uses

```text
epsilon_G- = epsilon_inv_G + r_G,
epsilon_O- = epsilon_inv_O + r_O,
```

before the same inverse/product propagation.  Thus both subspace uncertainty
and the measured zero-consistent leakage are retained.  No best-fit block,
Schur coefficient or sector-dependent carrier choice is allowed.

## Frozen spectral classifications

For either companion compute all right eigenvectors `Z`, their condition
number and the common Bauer--Fike bound

```text
epsilon_eig = kappa(Z) epsilon_C
              + 1000 eps m max(1,||C||),
```

where `m=50` or `30`.  A nonfinite/singular eigenvector matrix makes the cell
open.

For every eigenvalue `mu`, with `r=|mu|`:

```text
r < 1 - 100 epsilon_eig     CONTRACTING_RESOLVED
r > 1 + 100 epsilon_eig     EXPANDING_RESOLVED
|r-1| <= 10 epsilon_eig     UNIT_CONSISTENT
otherwise                   MODULUS_OPEN.
```

Record the complete complex eigenvalue list, labels, eigenvector condition,
ordered singular values, spectral radius and largest singular value.  Do not
infer long-time growth from a single update.

## Schedule comparisons

For each sector, variant and available companion:

- match even/odd complex eigenvalue multisets by the deterministic Hungarian
  minimum-total-distance assignment and record the maximum matched distance;
- compare ordered singular-value lists directly.

The comparison error is the sum of the two companion/eigenvalue errors plus
the floating floor.  Use

```text
distance <= 10 error    SCHEDULE_ROBUST
distance > 100 error    SCHEDULE_DEPENDENT
otherwise               SCHEDULE_OPEN.
```

There are `16` full-companion comparisons and, if every negative carrier is
invariant, `16` negative-companion comparisons.  Every available comparison
is recorded; a missing restricted companion is not silently removed but is
accounted for by the invariance outcome.

## Outcome ladder frozen before execution

Apply the first matching outcome:

1. `NEGATIVE_SHAPE_DYNAMICS_CONTROL_FAILED` if provenance, geometry, carrier
   or inherited `15+10` inertia contradicts the frozen inputs;
2. `NEGATIVE_SHAPE_CARRIER_OPEN` if the inherited carrier/gap is open;
3. `NEGATIVE_SHAPE_SUBSPACE_MIXED` if any of the 32 invariance residuals is
   `MIXING_RESOLVED`;
4. `NEGATIVE_SHAPE_INVARIANCE_OPEN` if none mixes but any invariance residual
   is open;
5. `NEGATIVE_SHAPE_COMPANION_SINGULAR` if any required `I+G` is
   singular-consistent;
6. `NEGATIVE_SHAPE_COMPANION_REGULARITY_OPEN` if any required regularity is
   open;
7. `NEGATIVE_SHAPE_SCHEDULE_DEPENDENT` if any available spectral or singular
   comparison is schedule-dependent;
8. `NEGATIVE_SHAPE_SCHEDULE_OPEN` if none is dependent but any is open;
9. `NEGATIVE_SHAPE_AUTONOMOUS_EXPANSION_RESOLVED` if the negative carrier is
   invariant and any restricted multiplier is expanding-resolved;
10. `NEGATIVE_SHAPE_AUTONOMOUS_MODULUS_OPEN` if the negative carrier is
    invariant but any restricted multiplier remains open;
11. `NEGATIVE_SHAPE_AUTONOMOUS_UNIT_CONSISTENT` if none is expanding-resolved
    or open and at least one restricted multiplier is unit-consistent;
12. `NEGATIVE_SHAPE_AUTONOMOUS_CONTRACTING_CENSUS` only if every restricted
    multiplier is contracting-resolved.

Full-sector expansion counts are reported under every outcome but do not
override the earlier invariance branches.

## Forbidden targets

The verifier must write `false` for loading or fitting any:

- desired expanding/contracting count;
- old tangent eigenvalue or eigenvector;
- scalar/vector/tensor harmonic spectrum;
- two-polarization target;
- continuum dispersion or speed;
- refinement, Planck or particle datum.

The only disclosed target is the already selected pair of negative-stiffness
sectors.

## Post-execution semantic correction

This section was added after the first outcome was visible.  The original
last branch was named `NONEXPANDING_CENSUS` when every entry was contracting
or unit-consistent.  That name overclaimed: `UNIT_CONSISTENT` means that the
uncertainty interval contains unit modulus, not that it excludes a slightly
expanding true value.  The corrected ladder separates `UNIT_CONSISTENT` from
an all-contracting census.  No matrix, error, threshold, eigenvalue label or
preceding outcome branch changed.
