# Raw-residual nonlinear defect correction: first frozen result

Date: 2026-08-13

Prior-art gate: `81b1aa1`

Preregistered protocol: `2d695c6`

Frozen implementation: `9d76d02`

Verifier: `reproducible/verify_gravity_600cell_dust_nonlinear_defect.py`

Machine-readable record:
`reproducible/gravity_600cell_dust_nonlinear_defect.json`

Targeted run: **10/10 implementation checks passed**.  The full suite was
deliberately not run, following the current instruction to test only the
active calculation.

## Result

**OPEN NUMERICALLY.**  No complete-action-validated nonlinear continuation
was found, but the preregistered kill condition was not met because the
frozen scan did not resolve all grid points.

The exact preregistered hit fractions are:

| denominator | hits |
|---|---:|
| signed direction/parity cases | 0/16 |
| direction/parity pairs | 0/8 |
| phase contrasts | 0/4 |

These zero hit counts are **DERIVED COMPUTATIONAL within the frozen scan**.
They are not an exclusion of the continuous sector because transverse
localization remained unresolved.

## What the defect iteration did

The verifier reused all 16 stored cases and all 80 stored grid starts from
the first nonlinear result.  It used only the committed quotient Hessian and
the preregistered damping sequence; no Jacobian, optimizer or favorable
restart was introduced.

For even parity:

- 8/40 grid points passed both binary64 transverse gates, exactly one at
  `t=0` in each signed case;
- 32/40 stopped because none of the eleven allowed damping factors strictly
  reduced the raw transverse residual;
- the final raw transverse norm over all grid points ranged from
  `6.485e-13` to `1.169e-8`;
- all eight `t=0` points were sent to the independent 100-decimal
  complete-action audit, and 0/8 passed.

For odd parity:

- 0/40 grid points passed both transverse gates;
- all 40 stopped because no allowed damping strictly reduced the raw
  transverse residual;
- the final raw transverse norm ranged from `1.632e-11` to `1.394e-8`;
- consequently there was no candidate to send to the action audit.

Thus all eight even cases and all eight odd cases have outcome
`NONLINEAR_CONTINUATION_NUMERICALLY_UNRESOLVED` under the frozen rules.

## Independent action audit of the eight even candidates

All eight candidates passed the branch, imaginary-contamination,
absolute-correction and binary/action-agreement gates.  Every candidate
failed the same three zero-consistency gates:

| gate quantity | observed range | independent error range |
|---|---:|---:|
| full equation norm | `7.066e-13`--`1.101e-12` | `2.640e-14`--`2.640e-14` |
| transverse correction norm | `6.173e-8`--`1.646e-7` | `1.706e-12`--`1.707e-12` |
| absolute collective scalar | `9.754e-17`--`9.790e-17` | `8.532e-19`--`8.532e-19` |

The binary/action equation discrepancy was only
`4.044e-13`--`4.713e-13`, below the frozen `1e-8` gate, so this is not an
implementation disagreement.  The apparently small binary64 residual is
not evidence of a stationary solution once compared with the independently
estimated errors.

## Status ledger

- **DERIVED COMPUTATIONAL:** the registered verifier passes 10/10 and exactly
  reproduces the frozen cases, starts, damping rule and action audit.
- **DERIVED COMPUTATIONAL:** no tested candidate is a validated stationary
  point; the registered hit fractions are 0/16, 0/8 and 0/4.
- **OPEN:** the nonlinear relative-phase sector at `eta=1e-4` and
  `|t|<=0.1`, because 72/80 frozen grid solves did not meet the transverse
  gates.
- **PATTERN-informed method:** the raw-residual acceptance metric was chosen
  after inspecting the preceding one-step improvement.  Preregistration
  prevents further adaptation from being counted as blind evidence, but
  does not restore blind provenance.
- **NOT TESTED:** amplitude scaling, the other 25 boundary directions, a
  second slab and the full 840-edge carrier.
- **NOT DERIVED:** a physical vacuum, a clock, a causal speed limit, Planck
  time or Planck mass.

## Framing judgement

The fixed-Hessian iteration is a numerical continuation probe, not a law of
evolution.  Its failure cannot refute all dynamics on the 600-cell, while a
success would have established only a local discrete stationary branch.
The present result therefore does not advance the physical claim.  It does
show that the earlier dramatic one-step residual reduction was insufficient:
under the frozen complete-action audit it produced no stationary solution.

The next defensible numerical question is why the explicit defect iteration
stagnates--binary evaluation floor, rotation of the soft subspace, or a real
absence of a root.  Any modified solver must be treated as a new,
post-result-informed protocol and cannot retroactively strengthen this test.
