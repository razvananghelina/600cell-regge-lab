# Corrected adversarial protocol: scale-invariant pseudo-longitudinal residual

Date: 2026-08-19

This target-disclosed correction follows the frozen contradictory audit.  It
does not replace or reinterpret that artifact.  The known midpoints are
approximately `2.963e-5` before applying `B^-1` and `5.483e-8` afterward.

## Frozen provenance

In addition to the complete input set of the preceding adversarial protocol,
require:

| input | SHA-256 |
|---|---|
| `verify_gravity_600cell_dust_action_york_direct_precision_adversarial.py` | `719d46bc519e152c06fd0bc064962532ae759303fcfd69759c89dd5cd7bc8352` |
| `gravity_600cell_dust_action_york_direct_precision_adversarial.json` | `e39203741513f128a208f22896abef53daa12db089ee7e43abf9c90643fc579b` |
| `gravity_600cell_dust_action_york_direct_precision_adversarial_protocol.md` | `a11ff07846ff62335f9b029af1b680007e68b3d165e8aedb114d71dd907ce0a9` |
| `gravity_600cell_dust_action_york_direct_precision_adversarial_result.md` | `b9dd475651a1a840d7bfc07e32a8dab72c32f591b9db5e2614b06e667384624e` |

Require the preserved contradictory outcome, all sixteen resolved span
residuals, all sixteen zero-consistent inverse residuals and all sixteen
augmented ranks 24.

## Construction

Repeat the independent centered-archive, normalized binary rigidity and
column-pivoted QR construction literally.  Do not import the primary direct
verifier or change the carrier.  Require ranks `470/354/4`, the `5+25` action
split, the `15+10` longitudinal split and positive-definite `B`.

## Dimensionless observables

For every cell compute

```text
rho_span = ||(I-Q_BL Q_BL*) A L||_2 / ||A L||_2,
rho_comm = ||(I-L L*) B^-1 A L||_2 / ||B^-1 A L||_2.
```

Abort the identity verdict if either denominator is below `1e-12`.  Record the
two exact conditioning inequalities

```text
r_span <= ||B||_2 r_comm,
r_comm <= ||B^-1||_2 r_span,
```

with an additive binary floor.  These inequalities are controls relating the
two distances; neither is fitted to the observed values.

For each dimensionless scalar use

```text
epsilon_relative = maximum displacement over four derivative variants
                 + 1000 eps_binary * 30 * kappa_relative,
```

where `kappa_relative` is the maximum of the longitudinal QR condition,
`cond(B)`, `cond(BL)`, `||A||/||AL||` and
`||B^-1 A||/||B^-1 A L||`.  These ratios are scale invariant.  No
dimensionful norm of `A`, `B` or `D` multiplies this floor.  Apply the same
`10/100` bands.

The augmented rank uses a relative singular threshold

```text
100 epsilon_relative * sigma_max([B L,A L]).
```

## Outcome hierarchy

1. `SCALE_INVARIANT_ADVERSARIAL_CONTROL_FAILED` if provenance, upstream,
   rank, carrier, denominator, conditioning or inequality controls fail.
2. `SCALE_INVARIANT_DIRECT_REFUTATION_CONFIRMED` if both relative residuals
   are `NONZERO_RESOLVED` and the augmented rank is greater than 15 in all
   sixteen cells.
3. `SCALE_INVARIANT_DIRECT_REFUTATION_REFUTED` only if both relative
   residuals are `ZERO_CONSISTENT` and the augmented rank is exactly 15 in at
   least one otherwise well-conditioned cell.
4. `SCALE_INVARIANT_DIRECT_REFUTATION_OPEN` otherwise.

This is still a DERIVED COMPUTATIONAL / STRUCTURAL result rather than a
symbolic or formal interval theorem.  No full suite is run.
