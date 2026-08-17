# Disclosed protocol: gravity--dust origin of the lapse stiffness

Date: 2026-08-17  
Prior-art commit: `475f534`.

This is **not blind preregistration**.  The total Schur scale and a rough
gravity--dust cancellation ratio were inspected after the full-rank result.
The purpose of this protocol is to freeze the algebraic decomposition,
uncertainty rule and claim boundary before writing its verifier.  A passing
result remains **PATTERN / DERIVED DECOMPOSITION**, not a new prediction.

## 1. Frozen input

Use only

```text
reproducible/gravity_600cell_dust_full_lapse_schur.json
SHA256 4a441ce6b328ffcbb1b673e1c932d411c6a8a00434107bc010e44537190a9349
```

and the committed dust constants/state contained in its frozen upstream input.
Require `18/18`, outcome `FULL_LAPSE_SCHUR_REGULAR`, both parity ranks 120,
zero nullity and zero open directions.

Do not load a continuum value, experimental value, desired tick, wave speed or
Planck scale.

## 2. Exact algebraic split

Recompute at 100 decimal digits

```text
M = (90/pi) (2pi-5 acos(1/3)) L0,
rho = 0.0102^2 exp(r),
h_dust = -(2pi M/120) sqrt(rho),
```

where `r` is the accepted state and both parities must agree.

Because the dust Hessian is diagonal only in the 120 weak pole coordinates,
it changes only the `D` block of the frozen Schur partition.  For every stored
sector midpoint define exactly

```text
S_gravity = S_total-h_dust I,
S(mu) = S_gravity+mu h_dust I,
mu in {0, 1/2, 1, 2}.
```

No additional Hessian evaluation or fitted coefficient is permitted.

## 3. Scalar-consistency rule

For each `n=5d` sector define

```text
alpha = Tr(S_total)/n,
delta_scalar = ||S_total-alpha I||_2,
epsilon = stored epsilon_global
          + sqrt(n^2) * stored maximum_entry_radius.
```

Classify it `SCALAR_CONSISTENT` only if

```text
delta_scalar < 10 epsilon.
```

The same deviation applies algebraically to `S_gravity` and every `S(mu)`.
Require the maximum difference among all 14 complex `alpha` values to be below
ten times the largest sector epsilon before assigning a common scalar.

## 4. Cancellation diagnostics

If all sectors are scalar-consistent, define from the common real midpoint

```text
g = alpha-h_dust,
mu_star = -g/h_dust = 1-alpha/h_dust,
cancellation_ratio = |alpha|/max(|g|,|h_dust|).
```

Report `g`, `h_dust`, `alpha`, `mu_star`, `|1-mu_star|`, and the signs and
magnitudes of the four frozen `S(mu)` scalars.

The following post-observed descriptive bins are frozen only to prevent later
verbal inflation:

```text
cancellation_ratio < 1e-5       NEAR_CANCELLATION_PATTERN
1e-5 <= ratio <= 1e-2           PARTIAL_CANCELLATION
ratio > 1e-2                    NO_STRONG_CANCELLATION
```

Even `NEAR_CANCELLATION_PATTERN` is not evidence that the mass was selected by
this condition.  The mass normalization predates the calculation but was not
derived here from minimization of `cancellation_ratio`.

## 5. Mechanical outcomes

1. `LAPSE_STIFFNESS_ORIGIN_CONTROL_FAILED` for provenance, parity, reality,
   dimension or uncertainty failure;
2. `LAPSE_STIFFNESS_NONSCALAR` if any sector or the cross-sector common-scalar
   gate fails;
3. `LAPSE_STIFFNESS_SCALAR_NEAR_CANCELLATION` if all scalar gates pass and the
   descriptive bin is `NEAR_CANCELLATION_PATTERN`;
4. `LAPSE_STIFFNESS_SCALAR_NO_NEAR_CANCELLATION` for either other bin.

All four are valid scientific outcomes.  The verifier exits successfully when
the outcome is assigned mechanically and all implementation controls pass.

## 6. Claim boundary

A passing near-cancellation result establishes:

- **DERIVED ALGEBRAIC:** the dust contribution is a scalar affine shift of the
  already-certified pole Schur operator;
- **DERIVED COMPUTATIONAL:** the numerical sizes of the two terms and residual;
- **PATTERN:** scalar consistency and the near-cancellation description.

It does not establish an exact scalar theorem, a Hamiltonian constraint, a
dust clock, gauge symmetry, refinement stability, a graviton, a limiting speed
or a fundamental unit of time.
