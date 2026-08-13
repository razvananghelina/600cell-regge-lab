# Local-Jacobian audit of nonlinear defect stagnation

Date: 2026-08-13

Prior-art gate: `53968b0`

Preregistered protocol: `e53dcaf`

Frozen implementation: `46186d6`

Registered verifier:
`reproducible/verify_gravity_600cell_dust_stagnation.py`

Machine-readable record:
`reproducible/gravity_600cell_dust_stagnation.json`

Targeted run: **13/13 implementation checks passed**.  The full suite was
deliberately not run under the current instruction to execute only tests for
the active calculation.

## Frozen aggregate result

```text
MIXED_OR_OTHER_STAGNATION
```

All `80/80` displaced local Jacobian inverse actions passed the
preregistered resolution gates.  Their mechanism counts are:

| mechanism | count |
|---|---:|
| `FIXED_MODEL_MISMATCH` | 54/80 |
| `RESOLVED_BUT_NOT_FIXED_MISMATCH` | 26/80 |
| `PRECISION_LIMITED_LOCAL_MODEL` | 0/80 |

The frozen dominant-mismatch gate required at least `60/80`.  The observed
`54/80` therefore does **not** receive the stronger label, even though it is
a majority.

This result is **DERIVED COMPUTATIONAL** for the frozen finite diagnostic and
has **PATTERN-informed method provenance**, because the diagnostic was
selected after the fixed-Hessian run stagnated.

## What was refuted

The hypothesis that the earlier iteration stopped at a binary64 Jacobian
floor is refuted for the particular inverse actions tested here:

- all 80 response vectors were stable under the two frozen Richardson
  extrapolations;
- relative response changes ranged from `1.74e-4` to `6.87e-3`;
- the estimated Jacobian-error action on the response ranged from
  `5.09e-7` to `4.41e-3` of the residual;
- every one of the 16,320 displaced Jacobian evaluations stayed finite and
  on the certified Lorentzian branch.

There is no large exchange of the four-dimensional soft subspace: every
principal cosine between the local and background soft subspaces lies between
`0.9999999999977` and numerical unity.  However, converting those cosines to
angles gives a maximum principal angle of about `2.14e-6` radians at
`|t|=0.1`, about `1.05e-6` at `|t|=0.05`, and only order `1e-8` at the centre.
With condition number of order `5e10`, this microscopic rotation can be
load-bearing.  Therefore the broad claim that "rotation caused the failure"
is not isolated as a distinct mechanism, but neither is small amplified
rotation refuted.

The local soft curvatures do not change sign.  The background and all local
quotient Jacobians are positive definite, and their four smallest
eigenvalues remain approximately `3.67e-8`--`5.61e-8` over the scan.

## Actual fixed-model failure

The matrix frozen at `t=0` was transported across the whole collective
interval.  Relative to that matrix, the local Jacobian changes in spectral
norm by approximately:

| `t` | relative change |
|---:|---:|
| `-0.10` | `5.127%` |
| `-0.05` | `2.532%` |
| `0` | `0.0017%`--`0.0021%` |
| `+0.05` | `2.469%` |
| `+0.10` | `4.877%` |

Both matrices remain positive definite, but their condition number is of
order `5e10`.  Consequently a small noncommuting change can make
`-H_Q^-1 F` increase `norm(F)^2/2`; positive definiteness of each matrix does
not make their nonsymmetric product a descent operator.  The change includes
both curvature variation within the soft sector and the small soft-subspace
rotation quantified above.

For even parity the old direction was robust non-descent at `31/40` states;
for odd parity it was non-descent at `28/40`.  Exactly one local-Jacobian
trial was allowed:

- any residual reduction occurred at `72/80` states;
- a reduction by more than one half occurred at `59/80` states;
- at noncentral `t`, that strong reduction occurred at `59/64` states;
- all `16` central states failed the strong-reduction gate, accounting for
  most of the aggregate shortfall.

The noncentral breakdown is a post-result diagnostic and is therefore
**PATTERN**, not a replacement for the frozen `54/80` outcome.

## Complete-action anchors

The deterministic audit reconstructed the full equation at 100 decimals for
all 16 signed cases at `t=-0.1,0,+0.1`: 48 anchors and 10,080 action
evaluations.  All 10,080 displaced geometries passed the branch gates,
maximum imaginary contamination was `9.151e-93`, and the maximum norm
between the binary analytic and complete-action rows was `1.015e-10`, below
the frozen `1e-8` ceiling.

The preregistered transverse classifications were:

| parity and location | result |
|---|---:|
| even, both endpoints | 16 resolved nonzero |
| even, centre | 8 unresolved |
| odd, endpoints and centre | 24 resolved nonzero |
| all zero-consistent anchors | 0 |

Thus the totals are `40/48 ACTION_RESOLVED_NONZERO`,
`8/48 ACTION_RESIDUAL_UNRESOLVED` and `0/48 ACTION_ZERO_CONSISTENT`.
No stored state is certified stationary.

The action-derived local correction norms are typically below `1e-5`, but
the binary/action bias, after inversion through the soft sector, is often
comparable to the correction.  A continuation which uses the binary
equation as its residual can therefore converge to the wrong numerical
root.  The next calculation must use the complete-action residual, even if a
binary local Jacobian is retained only as a preconditioner.

## Post-result scalar pattern

After the frozen result, the collective scalar was projected from the stored
complete-action rows at all 48 anchors.  It is positive in every case:

```text
8.21e-17 <= g <= 1.20e-16,
96.7 <= |g|/empirical_error <= 148.0.
```

This target was not part of the stagnation protocol.  The sign-definite
result is therefore **PATTERN**, not evidence of a derived nonlinear
obstruction.  Moreover the anchor states still have nonzero transverse
equations.  The scalar must be recomputed only after those equations are
solved with the complete-action residual.

The pattern is consistent with the known Regge mechanism in which a
linearized lapse freedom becomes a higher-order pseudo-constraint.  It is
not a novel mechanism and is not yet established on this carrier.

## Status ledger

- **DERIVED COMPUTATIONAL:** 13/13 targeted checks; 80/80 resolved local
  inverse actions; exact frozen mechanism and action-anchor counts.
- **DERIVED NEGATIVE:** a particular-response binary64 precision floor does
  not explain the old stagnation.
- **DERIVED:** there is no large soft-subspace exchange, but a smooth
  order-`1e-6` rotation away from the centre remains potentially material
  after amplification by the condition number.
- **DERIVED COMPUTATIONAL:** the frozen background Hessian is an inadequate
  global nonlinear solver on much of the collective interval.
- **PATTERN:** noncentral local-Newton robustness and the positive collective
  scalar at all 48 anchors.
- **OPEN:** whether complete-action transverse roots exist and whether their
  reduced collective scalar vanishes.
- **OPEN / NOT TESTED:** root nonexistence.  A failed solver is not a no-root
  certificate.
- **NOT TESTED:** amplitude scaling, the other 25 boundary directions, a
  second slab and the full 840-edge carrier.
- **NOT DERIVED:** a vacuum, physical time, inertia, mass, a causal speed
  limit or Planck units.

## Next decision

The next falsifiable object is the complete-action Lyapunov--Schmidt reduced
scalar.  For every frozen case and collective grid point, solve the 34
transverse equations using the 100-decimal action row as the residual and a
recomputed local Jacobian only as a preconditioner.  Then evaluate the
collective scalar and its error.  A sign-definite nonzero scalar for both
boundary signs after all transverse solves would be a genuine local
pseudo-constraint obstruction within the frozen scan.  A zero would produce
a nonlinear stationary candidate.  Either outcome requires a new
preregistered protocol.
