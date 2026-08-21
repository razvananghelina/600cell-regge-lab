# Protocol: combined envelope for constrained action/Hessian comparison

Date: 2026-08-21

This algebraic correction audit is frozen after the directional diagnostic in
commit `956ccf2` and before computing any combined comparison.

## 1. Exact question

The direct action ladder has an action-only uncertainty `e_action`; the stored
constrained response matrix has an independent certified entrywise uncertainty
`e_K`.  For a frozen coefficient vector `y`, test the same scalar equality
using the propagated combined bound rather than treating the finite-difference
Hessian as exact.

No action, Hessian, solve, class census or geometry is recomputed.

## 2. Frozen inputs

Require:

```text
reproducible/gravity_600cell_refined_h4_constrained_response.json
  f029260c9ee6e3b763293d237aae27e6ff7c1256eb8bc19c35725084ff385888
reproducible/gravity_600cell_refined_h4_constrained_directional_diagnostic.json
  35662f71e4debdbd64356c6e004d32f652719baf17b38843414bb25b96e21b58
docs/gravity/gravity_600cell_refined_h4_constrained_response_primary_first_result.md
  633a57f3d2b4a054cce20d08544d409dac8fdaf53c39bae72ab2e9fceb4e83eb
docs/gravity/gravity_600cell_refined_h4_constrained_directional_diagnostic_result.md
  0888a9b5caad440b1643d5f63992631a405f25ae13850888c3df0cc305a5ecb8
```

Require the primary `18/19 CONTROL_FAILED` artifact with one tentative class,
and the corrected `15/15 NONASYMPTOTIC` diagnostic with twelve complete
records.  Match records by exact schedule order and direction label.

## 3. Propagated comparison

For each direction read:

```text
q          = stored constrained quadratic response,
q_action   = final 180-digit Y_1 estimate,
e_action   = stored diagnostic envelope,
e_K        = primary entrywise response-matrix envelope.
```

The exact componentwise inequality

```text
|y^T Delta K y|
 <= sum_ab |y_a| |Delta K_ab| |y_b|
 <= ||y||_1^2 e_K
```

fixes, without a fitted multiplier,

```text
e_hessian=||y||_1^2 e_K,
e_total=e_action+e_hessian.                       (1)
```

The frozen coefficient one-norms are `1` for `first_basis_vector` and `11`
for `all_ones` and `alternating_signs`.  Require all twelve

```text
|q_action-q| <= e_total.                          (2)
```

Report every ratio to both `e_action` and `e_total`, not only the maximum.

## 4. Action convergence and controls

Independently of `q`, reconstruct from the stored diagnostic ladder

```text
|R0-R1|/|R1-R2|,
|R1-R2|/|R2-R3|,
|X0-X1|/|X1-X2|.
```

Require the first two ratios in `[8,32]` and the last in `[32,128]` for every
direction, bracketing the theoretical values `16,16,64`.

Controls:

1. Removing `e_hessian` must reject all twelve comparisons; this proves that
   the correction is load-bearing and not a redundant relabelling.
2. Replace `q` by `q_bad=q+1e-6*max(1,|q|)` and require the corruption to lie
   more than `10^6 e_total` away from `q_action` for every direction.
3. Require `e_action>0`, `e_K>0`, and the recorded matrix/class identities to
   remain byte-frozen; no envelope is recomputed from the observed mismatch.
4. No geometry, action, Hessian, solve, schedule classification, spectrum,
   continuum target or physical constant is evaluated.

## 5. Frozen outcomes

Use the first applicable outcome:

1. `REFINED_H4_COMBINED_ENVELOPE_CONTROL_FAILED` for provenance, record
   matching, convergence, positivity, action-only rejection, corruption or
   scope failure;
2. `REFINED_H4_COMBINED_ENVELOPE_MISMATCH` if any of the twelve comparisons
   fails (2);
3. `REFINED_H4_COMBINED_ENVELOPE_CORROBORATED` otherwise.

Outcome 3 licenses a separately preregistered primary-verifier correction:
replace only its directional scalar check by the same high-order action ladder
and bound (1), then rerun the complete primary verifier twice.  It does not
accept the tentative schedule class on its own.

## 6. Execution

Register before first execution, run twice with a byte-identical artifact,
and perform only the static registry audit.  Do not run the full suite.

