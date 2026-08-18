# Preregistered blind protocol: complete two-slab tangent cocycle

Date: 2026-08-18

Prior-art gate commit: `d51fbbf`.

This is a blind operator census.  No second-slab full tangent, product
spectrum, amplification count, `S^3` harmonic target or desired physical
classification has been inspected before this protocol commit.

## Frozen inputs

| input | SHA-256 |
|---|---|
| `gravity_600cell_dust_full_boundary_tangent.json` | `4da8bcd2890a54bc9d3b60c6195df2933ea56194d942ab0285b51599ba287bd5` |
| `gravity_600cell_dust_full_boundary_tangent.npz` | `816c605da2a655442bbadce7a23965f0822f99e7bdc1d0a4a27af548de85446b` |
| `verify_gravity_600cell_dust_full_boundary_tangent.py` | `c8662bb0835865aac6696fc3f474ed668fed3fe393b9c32a59e709a984c35571` |
| `gravity_600cell_dust_homothetic_canonical_lapse.json` | `4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9` |
| `gravity_600cell_dust_second_tick_local_correction.json` | `936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70` |
| `verify_gravity_600cell_dust_second_tick_local_correction.py` | `cef59fa0bc3a1c8fa3be0193234371b7dda303a0ec72683ddcdd88bcb40f3725` |
| `verify_gravity_global_regge_orbits.py` | `ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf` |
| `verify_gravity_600cell_dust_full_anisotropic_legendre_rank.py` | `834b97c85e386def853b6308e65e831c52d62d7cbcc4b23118602120d6c676e5` |

The first tangent archive contains 224 midpoint/radius arrays.  Its binary64
midpoints are re-enclosed before multiplication by adding one half-ULP in
each real and imaginary component to the stored Flint radius.  This prevents
the new product ball from pretending that binary serialization was exact.

Both schedule parities and the four already frozen derivative variants are
mandatory.  Arithmetic uses 100 decimal digits and 80-decimal Flint balls.

## Exact second-slab geometry

For each parity use its own committed values, without averaging:

```text
a_old = a_1,
a_new = a_2,
r     = r_2,

q_old      = exp(2 a_1) L_0^2,
q_new      = exp(2 a_2) L_0^2,
q_diagonal = exp(a_1+a_2) L_0^2-rho_2,
q_pole     = -rho_2.
```

The dust Hessian uses the original conserved total `MASS`.  Require the
complete committed thirty-component first/second seam, all base and displaced
Lorentzian branch gates, reality, free-action carrier, and the literal
boundary layer map before constructing a tangent.

## Complete block reconstruction

Derive the same seven deterministic minimal sector bases from the literal
group action.  For every parity, sector and derivative variant:

1. project the complete second-slab Hessian;
2. require the full `35d+30d` pre-Legendre determinant ball to exclude zero;
3. solve the canonical response and form `T_2`;
4. verify the complex-sector identity `T_2* Omega T_2=Omega` with ball
   arithmetic;
5. re-enclose the committed matching `T_1` block and form the rigorous ball
   product `C_21=T_2T_1`;
6. verify `C_21* Omega C_21=Omega` directly from the product ball.

The sector weights are fixed as `d`, giving

```text
sum_sector d * (60d) = 1440.
```

No minimal-block count may be reported as a full-carrier count without this
weight.

## Blind spectral census

For `T_2` and `C_21` separately, and for all four derivative variants, report:

- ordered singular values;
- eigenvalues and reciprocal-conjugate matching defect;
- determinant modulus/log modulus;
- tangent condition number and eigenvector condition number;
- symplectic defect and its propagated error;
- number of eigenvalues labelled unit-consistent, resolved off-unit or open;
- resolved contracting and expanding counts, restored with representation
  multiplicity;
- minimum and maximum resolved moduli where they exist.

Reuse the already corrected SVD calibration.  If `sigma_max` is large, the
absolute backward-error floor is proportional to
`eps_machine*sigma_max`; the extreme reciprocal-product floor is therefore
proportional to `eps_machine*sigma_max^2`, not machine epsilon alone.

An eigenvalue is

```text
UNIT_CONSISTENT       ||lambda|-1| < 10 epsilon_lambda,
RESOLVED_OFF_UNIT     ||lambda|-1| > 100 epsilon_lambda,
NUMERICALLY_OPEN      otherwise,
```

where `epsilon_lambda` includes all derivative variants, the full product
ball radius, non-normal eigenvector conditioning and binary64 floor.  Open
eigenvalues are never silently allocated to either branch.

## Schedule comparison

The even/odd minimal matrices need not be entrywise equal in their internal
orbit bases.  Compare target-independent spectral invariants sector by
sector.

The primary schedule diagnostic is the maximum distance between the ordered
singular spectra.  Its uncertainty is the sum of both maps' singular-value
uncertainties plus serialization floor:

```text
SCHEDULE_ROBUST       distance <= 10 epsilon,
SCHEDULE_DEPENDENT    distance > 100 epsilon,
SCHEDULE_OPEN         otherwise.
```

Optimally matched eigenvalue distance is reported with the analogous labels
but is secondary because the maps may be highly non-normal.  It cannot veto a
singular-spectrum result unless it is itself resolved dependent under its
larger calibrated uncertainty.

## Frozen outcome hierarchy

1. `TWO_STEP_FULL_TANGENT_CONTROL_FAILED` for any provenance, seam, branch,
   carrier, sector-exhaustion, determinant, symplectic or archive failure.
2. `TWO_STEP_FULL_TANGENT_SCHEDULE_DEPENDENT` if any `T_2` or `C_21` sector is
   resolved dependent by the primary singular-spectrum test, or by a resolved
   eigenvalue test not contradicted by its conditioning gate.
3. `TWO_STEP_FULL_TANGENT_SCHEDULE_OPEN` if no sector is dependent but at
   least one primary singular comparison is open.
4. `TWO_STEP_FULL_TANGENT_COCYCLE_CERTIFIED` only if all fourteen primary
   comparisons (`7` sectors times `T_2,C_21`) are robust and every structural
   control passes.

The outcome is independent of how many eigenvalues expand or contract.

## Explicit exclusions

- no comparison with continuum `S^3` scalar/vector/tensor harmonics;
- no assignment of gauge, graviton or matter labels;
- no physical energy norm, Lyapunov exponent, frequency or limiting speed;
- no third tick, refinement or nonlinear anisotropic solve;
- no full-suite run.
