# Preregistration: local-Jacobian audit of nonlinear defect stagnation

Date: 2026-08-13

Prior-art gate: `53968b0`

Frozen nonlinear defect result: `8f651cd`

Status: **frozen before evaluating any displaced local Jacobian, trying any
local-Jacobian step, or evaluating a new complete-action anchor**.

This is a post-result-informed diagnostic.  It cannot be counted as a blind
discovery test and it is not a solver for a root.

## 1. Frozen carrier, cases and states

Load without alteration:

- the two schedule parities;
- the four stored Helmert phase contrasts and both signs;
- `eta=1e-4` and all stored final-boundary logarithmic vectors;
- the five grid values `t=-0.10,-0.05,0,+0.05,+0.10`;
- the final `z` written by the frozen defect result for every one of the 80
  grid states;
- the exact complement `Q` and parity-specific committed background quotient
  Hessian `H_Q`.

Use the unchanged complete analytic logarithmic equation

```text
E = u*(partial S_total/partial u)/24,
F(z;t,b) = Q^T E(log u_path(t)+Qz,b).
```

No restart, amplitude change, boundary direction, lapse value or favorable
state may be selected.

## 2. Displaced local Jacobian

At every stored state compute centered differences in each of the 34 `z`
coordinates at the frozen logarithmic steps

```text
h1=5e-4, h2=2.5e-4, h3=1.25e-4.
```

Every displaced evaluation must remain finite and pass the certified
Lorentzian branch gates.  Define

```text
J1=J(h1), J2=J(h2), J3=J(h3),
R12=(4*J2-J1)/3,
R23=(4*J3-J2)/3,
J6=(16*R23-R12)/15,
dJ=J6-R23.
```

Because `F` is the gradient of the restricted action in logarithmic
coordinates, its exact Jacobian is symmetric.  Record

```text
symmetry_defect = norm_2(J6-J6^T)/2,
J = (J6+J6^T)/2,
error = norm_2(dJ).
```

Do not use the full-matrix error norm as a smallest-singular-value gate: it
can be dominated by truncation error in the thirty hard modes while the
particular weak response `J^-1 F` is stable.  Resolution of the inverse action
is therefore tested directly in section 4.

For diagnostics, also report all singular values, the four-dimensional
small-absolute-eigenvalue subspace of `J`, and its four principal cosines
with the corresponding four-dimensional soft subspace of `H_Q`.  This
subspace comparison is target-free and is not itself a pass gate.

## 3. Frozen-model forcing and descent

At every state form

```text
p_H = H_Q^-1 F,
eta_H = norm(F-J*p_H)/max(norm(F),1e-30),
mu_H = F^T J p_H/max(norm(F)^2,1e-60),
delta_mu_H = norm(dJ*p_H)/max(norm(F),1e-30).
```

Here the actual update direction is `-p_H`; positive `mu_H` is infinitesimal
descent for `norm(F)^2/2`, while negative `mu_H` is ascent.  Assign:

- `ROBUST_FIXED_NONDESCENT` if
  `mu_H+10*delta_mu_H < 0`;
- `ROBUST_FIXED_DESCENT` if
  `mu_H-10*delta_mu_H > 0`;
- `FIXED_DIRECTION_UNRESOLVED` otherwise.

These labels are used only when the local Jacobian gate passes.  Record the
full counts rather than selecting an illustrative state.

## 4. Exactly one local-Jacobian trial

Form also

```text
J23=(R23+R23^T)/2,
p_J=J^-1 F,
p_23=J23^-1 F,
step_change=norm(p_J-p_23)/max(norm(p_J),1e-30),
step_model_error=norm(dJ*p_J)/max(norm(F),1e-30).
```

A local inverse action is `JACOBIAN_STEP_RESOLVED` only if:

- every displaced branch passes;
- `symmetry_defect <= 10*max(error,1e-30)`;
- both symmetrized matrices have smallest singular value above the frozen
  absolute numerical floor `1e-12`;
- `step_change<=0.1` and `step_model_error<=0.1`.

Otherwise label it `LOCAL_JACOBIAN_UNRESOLVED` and do not use the inverse as
evidence.  This is a certificate for the particular response, not a claim
that every entry or eigenvalue of the full matrix is resolved.

When this gate passes, try the same damping factors

```text
1,1/2,1/4,...,1/1024
```

in that order and retain only the first finite, branch-valid trial which
strictly reduces `norm(F)`.  Stop after this one accepted step.  Do not
iterate, localize the collective scalar or validate a root.

Record the accepted damping and reduction factor.  A state receives
`LOCAL_JACOBIAN_DESCENT` only if a trial is accepted and the reduction factor
is below `0.5`; otherwise it receives `NO_STRONG_LOCAL_JACOBIAN_DESCENT`.

The diagnostic mechanism label is:

- `FIXED_MODEL_MISMATCH` if the Jacobian is resolved, the old direction is
  `ROBUST_FIXED_NONDESCENT`, and the new direction has
  `LOCAL_JACOBIAN_DESCENT`;
- `RESOLVED_BUT_NOT_FIXED_MISMATCH` for any other resolved Jacobian;
- `PRECISION_LIMITED_LOCAL_MODEL` for an unresolved Jacobian.

No label means that a root exists or does not exist.

## 5. Independent complete-action anchors

Independently reconstruct the full 35-component equation from the complete
action at 100 decimal digits using logarithmic steps

```text
2e-5,1e-5,5e-6
```

and the already certified two-stage Richardson extrapolation.  The anchor
set is frozen as all 16 signed direction/parity cases at exactly

```text
t=-0.10, 0, +0.10,
```

for 48 anchors total.  This selection covers both scan endpoints and the
centre without looking at residual size.  All 210 displaced geometries per
anchor must pass the branch gates; maximum imaginary contamination must be
below `1e-80`.

For `F_action=Q^T E_action` and its projected empirical error `epsilon_F`,
assign:

- `ACTION_ZERO_CONSISTENT` if
  `norm(F_action)<=10*max(epsilon_F,1e-30)`;
- `ACTION_RESOLVED_NONZERO` if
  `norm(F_action)>100*max(epsilon_F,1e-30)`;
- `ACTION_RESIDUAL_UNRESOLVED` otherwise.

Also record `norm(F_binary-F_action)` and its ratio to `epsilon_F`; do not
relax the analytic/action agreement ceiling `1e-8`.

## 6. Aggregate outcomes

Report exact multisets and fractions over all 80 local states and all 48
action anchors.  Use one aggregate outcome:

- `FIXED_MODEL_MISMATCH_DOMINANT` if at least 60/80 states have resolved
  Jacobians and at least 3/4 of those resolved states receive
  `FIXED_MODEL_MISMATCH`;
- `LOCAL_MODEL_PRECISION_LIMITED` if fewer than 60/80 Jacobians resolve;
- `MIXED_OR_OTHER_STAGNATION` otherwise.

The `60/80` coverage and `3/4` fraction are frozen before the matrices are
computed.  They are descriptive robustness thresholds, not physical
constants.

## 7. Claim boundary

- **DERIVED COMPUTATIONAL:** only the recorded finite matrices, residuals,
  branch checks, action anchors and frozen classification counts.
- **PATTERN-informed:** the diagnostic method, because it was chosen after
  seeing the fixed-Hessian stagnation.
- **OPEN:** a stationary root, a continuous nonlinear branch and absence of
  either.
- **NOT TESTED:** a second slab, the full 840-edge carrier, amplitude scaling
  or the other 25 boundary directions.
- **NOT DERIVED:** physical time, a vacuum, inertia, mass, a causal speed
  limit or Planck units.

If the local model is precision-limited, the next legitimate step is an
arbitrary-precision Jacobian or interval formulation.  If fixed-model
mismatch is dominant, a separately preregistered local-Jacobian continuation
may be attempted.  Neither successor may reuse this diagnostic as blind
evidence.
