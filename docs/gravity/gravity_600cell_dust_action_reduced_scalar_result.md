# Complete-action reduced-scalar scan: numerical resolution boundary

Date: 2026-08-14

Prior-art update: `8be130e`

Frozen stagnation result: `1d66278`

Preregistered protocol: `17f9560`

Frozen implementation: `d39a51a`

Operational interruption-safe checkpoint amendment: `27985ba`

Registered verifier:
`reproducible/verify_gravity_600cell_dust_action_reduced_scalar.py`

Machine-readable record:
`reproducible/gravity_600cell_dust_action_reduced_scalar.json`

Targeted run: **12/12 implementation checks passed**.  The full suite was
not run, following the standing instruction to execute only tests for the
active calculation.

## Headline

The frozen complete-action Lyapunov--Schmidt scan is

```text
ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED.
```

None of the 80 frozen grid states passed the independent transverse
validation.  Consequently no reduced collective scalar was admissible for
classification, no sign bracket was licensed, and there is no stationary
hit or physical obstruction in this result.

The calculation nevertheless isolates a sharp numerical mechanism.  The
solver and validation extrapolants agree to about `3.5e-16--5.7e-16`, but
the solver-window Richardson error proxy is almost exactly 81 times the
validation-window proxy.  It therefore declares nonzero residuals to be
zero-consistent and creates a line-search dead zone.  This falsifies the
frozen solver certification; it does not falsify a stationary continuation.

## Frozen aggregate result

The attempt and hit counts are

| quantity | result |
|---|---:|
| frozen grid states | 80 |
| mechanically forced bisection points | 0 |
| total attempts | 80 |
| signed-case hits | 0/16 |
| direction/parity-pair hits | 0/8 |
| phase-contrast hits | 0/4 |
| independently transverse-validated states | 0/80 |
| admissible scalar zero/nonzero/unresolved labels | 0/0/0 |

All sixteen signed cases therefore receive
`ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED`.

The solver outcomes before independent validation are

| parity | solver-zero | no robust action descent | accepted-iteration census |
|---|---:|---:|---|
| even | 40/40 | 0/40 | 26 once, 14 twice |
| odd | 23/40 | 17/40 | 7 zero, 11 once, 22 twice |

Here "zero" in the final column means no accepted correction, not a zero
equation.  The seven such odd states all ended with
`NO_ROBUST_ACTION_DESCENT`.

## Independent validation refutes all 63 solver-zero labels

Every one of the 63 solver-zero candidates passes the branch, imaginary-part
and frozen absolute binary/action agreement gates.  Every one fails the
independent complete-action transverse-zero gate.

| parity | candidates | fail `F` zero | fail preconditioned `p` zero |
|---|---:|---:|---:|
| even | 40 | 40/40 | 27/40 |
| odd | 23 | 23/23 | 19/23 |

The remaining gates pass for all candidates:

- propagated preconditioner error is below `1e-5`;
- all 210 validation geometries remain on the certified Lorentzian branch;
- imaginary contamination remains below `1e-80`;
- the binary/action norm difference remains below the preregistered absolute
  ceiling `1e-8`.

The decisive residual statistics are

| quantity | even min / median / max | odd min / median / max |
|---|---|---|
| validation `norm(F)` | `4.31e-13 / 7.06e-13 / 1.12e-12` | `4.03e-13 / 7.21e-13 / 1.06e-12` |
| validation error proxy | `2.27e-14 / 2.64e-14 / 3.07e-14` | `2.27e-14 / 2.45e-14 / 3.07e-14` |
| validation `norm(F)/error` | `15.83 / 26.42 / 40.95` | `13.79 / 26.94 / 43.40` |
| validation `norm(p)/error` | `0.753 / 49.92 / 748.20` | `4.55 / 154.15 / 643.82` |

The solver had instead reported

| quantity | even min / median / max | odd min / median / max |
|---|---|---|
| solver `norm(F)/error` | `0.196 / 0.326 / 0.506` | `0.170 / 0.333 / 0.536` |
| solver error proxy | `1.84e-12 / 2.14e-12 / 2.48e-12` | `1.84e-12 / 1.98e-12 / 2.48e-12` |

Thus the residual value is stable but its conservative proxy changes scale:

```text
solver error / validation error
  even: 80.9783 / 80.9804 / 80.9822,
  odd:  80.9783 / 80.9813 / 80.9822.
```

This is the expected factor `3^4=81` for the fourth-order Richardson
difference when the two frozen step windows differ by a factor of three.
The sixth-order extrapolated rows themselves differ by only
`3.47e-16--5.72e-16`, and their transverse norms agree within about
`0.025%`.  The nonzero validation residual is therefore not a disagreement
between the two action calculations.

**DERIVED COMPUTATIONAL:** the larger frozen window overestimates the error
proxy by about 81 relative to the smaller window and terminates too early.
This is a numerical certification failure, not a physical zero.

## The line-search dead zone

The 17 odd `NO_ROBUST_ACTION_DESCENT` states split into:

- seven states with zero accepted iterations;
- ten states with one accepted iteration.

For the latter ten, the current solver-window ratio was only
`norm(F)/error=0.228--0.617`.  The regular acceptance inequality was

```text
norm(F_trial)+10 error_trial
  < norm(F_current)-10 error_current.
```

Its right-hand side is negative whenever `norm(F_current)<=10 error_current`,
while its left-hand side is nonnegative.  Regular descent was therefore
mathematically impossible for all ten, independently of the proposed
direction.  They could be accepted only through the provisional-zero path,
whose preconditioned-zero gate failed on every damping.

The same impossibility already affects one of the seven zero-iteration
states.  In total, `11/17` no-descent labels begin inside this dead zone.

The other evidence does not support calling the preconditioner direction
uniformly bad:

- for the seven zero-iteration failures, the best raw trial residual was only
  `1.23%--3.68%` of the current residual;
- for the ten one-iteration failures, nine of ten best raw trials reduced the
  residual, with ratios `0.445--0.807`; one had ratio `1.039`;
- all were rejected because the frozen ten-error bands overlapped or because
  the provisional preconditioned zero was not certified.

**DERIVED NUMERICAL NEGATIVE:** `NO_ROBUST_ACTION_DESCENT` here cannot be
read as evidence that no root exists, or even as evidence that every raw
preconditioned direction is non-descent.

The exact 17 identities `(direction, boundary sign, t)` are

```text
zero accepted iterations:
  (1,-1,0), (1,+1,0), (2,-1,0), (2,+1,0),
  (3,+1,0), (4,-1,0), (4,+1,0)

one accepted iteration:
  (1,-1,+0.05), (1,+1,-0.10), (1,+1,+0.05),
  (2,-1,-0.05), (2,+1,-0.05), (2,+1,+0.10),
  (3,+1,-0.10), (4,-1,+0.10),
  (4,+1,-0.10), (4,+1,-0.05)
```

## Binary preconditioner boundary

The preregistered absolute gate `norm(E_binary-E_action)<1e-8` passes easily,
but it is too weak near a `1e-12` transverse residual.  On the 63 validation
states the absolute differences are

```text
even: 3.32e-13--5.42e-13,
odd:  2.54e-11--9.14e-11.
```

Relative to the complete-action transverse residual, those full-row norm
differences are `0.349--1.046` on the even parity and `29.1--216.2` on the
odd parity.  The comparison mixes a full-row difference with a transverse
norm, so it is a diagnostic rather than a relative-error theorem.  It is
nevertheless enough to show that the absolute `1e-8` gate cannot certify the
binary preconditioner at the scale needed for a root.

**DERIVED scope correction:** binary/action agreement at a coarse absolute
ceiling is not evidence that the binary Jacobian is an accurate local model
at the final residual scale.

## The positive scalar is inadmissible

All 63 raw validation scalars are positive:

```text
8.33e-17 <= g <= 1.03e-16,
102.6 <= abs(g)/error <= 114.7.
```

But every transverse validation failed first.  Under the frozen protocol all
63 labels are therefore `SCALAR_NOT_CLASSIFIED`.  The common positive sign is
at most a **PATTERN** consistent with the previously inspected anchors.  It
is not evidence for a pseudo-constraint and cannot be used to open a sign
bracket or claim a no-go.

## Hostile framing audit

1. `0/16` hits is not a stationary-point negative because no signed case was
   numerically resolved.
2. `0/80` transverse validations does not show that transverse roots do not
   exist.  It shows that the frozen stopping and validation windows disagree
   on whether the residual is negligible.
3. The exact factor near 81 is not a mysterious physical hierarchy.  It is
   the deterministic step-scaling of the Richardson error proxy.
4. The raw positive `g` values cannot be promoted after seeing them; all are
   inadmissible under the preregistered gate and retain PATTERN-informed
   provenance.
5. The 17 no-descent labels mainly diagnose the interaction between a large
   conservative error band and the acceptance rule.  They do not establish
   a bad physical direction or root nonexistence.

## Status ledger

- **DERIVED COMPUTATIONAL:** the registered targeted verifier passes
  `12/12` checks and writes all 80 frozen states.
- **DERIVED COMPUTATIONAL:** global outcome
  `ACTION_REDUCED_SCAN_NUMERICALLY_UNRESOLVED`.
- **DERIVED COMPUTATIONAL:** even solver outcomes `40/40` apparent zeros;
  odd outcomes `23/40` apparent zeros and `17/40` no robust descent.
- **DERIVED COMPUTATIONAL NEGATIVE:** independent validation rejects all
  `63/63` apparent transverse zeros.
- **DERIVED NUMERICAL MECHANISM:** the solver error proxy is about 81 times
  the validation proxy while the extrapolated action rows agree near
  `1e-16`.
- **DERIVED NUMERICAL MECHANISM:** at least `11/17` no-descent states lie in
  an algebraic dead zone of the frozen robust-descent inequality.
- **DERIVED SCOPE LIMIT:** the absolute binary/action agreement gate is not
  scale-sensitive enough near the candidate roots.
- **PATTERN:** every raw scalar available after the solver is positive.
- **OPEN NUMERICALLY:** existence of complete-action transverse roots and
  the sign or zero of the reduced scalar at them.
- **NOT DERIVED:** a stationary vacuum, physical time, inertia, mass, a
  causal speed limit, Planck units, or a particle-mass operator.

## Next falsifiable calculation

A new protocol may not relax the present result retroactively.  It should
freeze, before execution:

1. the current smaller validation window as the next solver window;
2. a third, disjoint and still smaller validation window;
3. a scale-sensitive certification of any binary preconditioner, based on
   its action on the proposed correction rather than the old absolute
   `1e-8` row difference;
4. a deterministic fallback to a complete-action Jacobian or trust-region
   step for states where that certification fails;
5. the same 80 states, signs, directions, grid and hit accounting.

That continuation would test a solver repair suggested by a fully disclosed
negative result.  It would remain PATTERN-informed in method provenance and
would still say nothing about masses until a stationary background and a
separately derived fermionic operator exist.
