# Directional diagnostic result: action converges past the stored Hessian value

Date: 2026-08-21

## Provenance and reproducibility

| stage | commit |
|---|---|
| diagnostic protocol | `0baeeee` |
| verifier registered | `42b3b88` |
| textual parser failure preserved | `dd7e05a` |
| parser correction preregistered | `1580829` |
| parser correction implemented | `678fb55` |
| default-precision failure preserved | `136da68` |
| precision correction preregistered/amended | `669917d`, `b769c4f` |
| precision correction implemented | `62abdd6` |

The corrected diagnostic passed `15/15` twice and produced the byte-identical
artifact

```text
reproducible/gravity_600cell_refined_h4_constrained_directional_diagnostic.json
SHA-256 35662f71e4debdbd64356c6e004d32f652719baf17b38843414bb25b96e21b58.
```

No Hessian, full suite or deferred nonlinear census was run.

## Frozen outcome

```text
REFINED_H4_DIRECTIONAL_DIAGNOSTIC_NONASYMPTOTIC
```

Only `4/12` Richardson error sequences remained monotone relative to the
stored primary quadratic value through every level, and `0/12` final
eighth-order action estimates lay inside the action-only envelopes around
that value.  This is the formal diagnostic result and is not overwritten.

## What the ladder itself shows

The high-precision reconstruction reproduces the original failed Richardson
values to `3.40e-59`.  For all twelve directions, independently of the stored
quadratic target,

```text
|R0-R1|/|R1-R2| = 16,
|R1-R2|/|R2-R3| = 16,
|X0-X1|/|X1-X2| = 64
```

to the displayed precision.  These are exactly the theoretical fourth- and
sixth-order factors.  The action ladder is therefore internally asymptotic;
the `NONASYMPTOTIC` label arose because its frozen gate measured convergence
*toward the finite-difference Hessian value*, treating that value as exact.

The final action/Hessian discrepancies are between approximately `3.85e-29`
and `8.78e-28`.  They are far above the action-only extrapolation envelopes,
but below the already frozen response-matrix uncertainty.  Using the rigorous
entrywise bound

```text
|y^T Delta K y| <= ||y||_1^2 e_K,
```

the maximum observed discrepancy is only about `6.64e-4` of the applicable
quadratic Hessian bound.  This last comparison was made after the diagnostic
outcome and is therefore **POST-DIAGNOSTIC STRUCTURAL**, not yet an accepted
replacement control.

## Interpretation

- **DERIVED COMPUTATIONAL:** the direct action differences have the exact
  expected step-halving factors through eighth order.
- **DERIVED NEGATIVE:** the stored finite-difference Hessian value cannot be
  treated as exact at the action extrapolation scale.
- **OPEN:** whether a separately frozen combined action-plus-Hessian envelope
  accepts all twelve comparisons.
- **NOT A RESULT:** the tentative single schedule class remains unaccepted
  while the primary verifier outcome is `CONTROL_FAILED`.

The next admissible step is a target-free algebraic audit of the combined
envelope `e_action+||y||_1^2 e_K`, with estimate-difference convergence gates
and a corruption control.  Only after that audit may the primary directional
criterion be corrected and the complete primary verifier rerun.

