# Preregistration: target-independent stationary momentum envelope

Date: 2026-08-16

Prior-art gate: `dedcbc6`.

Accepted geometry result: `46a7361`.

Status: **frozen before any stationary-curve evaluation**.

## 1. Target firewall

The verifier must not parse any momentum target or failed-solver artifact.
It may verify the SHA-256 of the accepted first-tick artifact as provenance,
then use only these already committed geometry constants:

```text
a1 = -3.11605957669450169173470644419863944122165192557277135128791e-6,
r1 = -3.55925313517063343725030533963917396571974345422547402551491e-6.
```

No desired pre-momentum may appear in the code, protocol or output.  The
comparison is a separate later mission.

## 2. Fixed curve and observables

With lower log fixed to `a1`, define for `(b,r)`

```text
q_old       = exp(2*a1)*L0^2,
q_new       = exp(2*b)*L0^2,
pole        = exp(r)*rho0,
diagonal    = exp(a1+b)*L0^2-exp(r)*rho0,
G(b,r)      = mean(five complete pole equations),
P(b,r)      = mean(thirty canonical pre-momenta).
```

The curve is `G=0`.  At every point compute the calibrated Jacobian

```text
J_GP = [[G_b,G_r],[P_b,P_r]],
D    = det(J_GP)=G_b*P_r-G_r*P_b.
```

On a regular graph `b(r)`, `D/G_b` is `dP/dr`; a resolved sign change of `D`
therefore brackets a momentum fold.

## 3. Base root without a target

At fixed `r=r1`, set

```text
b_center = 2*a1,
half_width = abs(a1).
```

Evaluate `G` at `b_center +/- half_width`; if their signs agree, double the
half-width, at most twelve times.  Accept the first Lorentzian resolved sign
bracket and perform exactly 100 bisections.  Require final width `<1e-30` and
`|G|<1e-25` at the midpoint evaluation.

No alternative center or bracket is allowed.

## 4. Frozen decreasing-lapse domain

Starting from the base root, solve `G=0` at the 32 further grid nodes

```text
r_j = r1-j/4,  j=1,...,32,
```

then at the four near-null sentinels

```text
r = r1-12, r1-16, r1-24, r1-32.
```

At each node use the preceding accepted `b` unchanged as the sole seed.  Apply
a scalar Newton solve with at most six accepted iterations and damping
`1,1/2,...,2^-10`, accepting the first branch-valid trial with

```text
|G_trial| <= (1-alpha/4)*|G_current|.
```

Stop successfully at `|G|<1e-25`.  No adaptive `r`, alternate seed, target,
optimizer or restart is permitted.

## 5. Derivative calibration

For both the scalar Newton derivative and `J_GP`, use central differences with

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Use the existing entrywise factor-10 cross-calibration.  For a scalar
derivative set its error to the sum of the three absolute cross-scale
differences plus `1e-60` and require `|G_b|>100*error` before Newton.

For `D`, compute all four determinant estimates and set

```text
epsilon_D = |Dop1-Dop2|+|Dval1-Dval2|+|Dop1-Dval1|+1e-60.
```

A determinant sign is resolved only if `|Dop1|>100*epsilon_D`.

Every perturbed evaluation must retain the complete Lorentzian/complex-angle
branch.

## 6. Full node and parity gates

At every curve node require

```text
max abs(30 diagonal equations) < 1e-60,
max abs(5 pole equations)      < 1e-25,
within-type spreads            < 1e-60,
resolved nonzero G_b,
resolved determinant sign unless the node lies inside a refined fold bracket.
```

Solve the curve with the even schedule.  Substitute the same `(b,r)` into the
odd schedule and require

```text
|G_even-G_odd| < 1e-24,
|P_even-P_odd| < 1e-22,
max abs(local_even-local_odd) < 1e-24.
```

At every refined fold midpoint also compute the complete odd derivative
matrix and require agreement of all entries within their combined calibrated
errors.

## 7. Fold enumeration and refinement

Enumerate every adjacent main-grid pair (`j=0,...,31`) with opposite resolved
signs of `D`.  Sentinel gaps do not define refinement brackets.  Record the
exact number before any refinement.  More than four brackets is an explicit
enumeration failure.

For each bracket perform exactly 24 bisections in `r`.  At every midpoint solve
`G=0` using the linear interpolation of the two bracket endpoint `b` values as
the sole seed, then compute the resolved determinant sign.  Retain the half
with opposite endpoint signs.  Require final width `<2e-8` and resolved
opposite signs at both ends.

Evaluate the midpoint curve state once more and record its momentum as the
fold momentum estimator.  Do not compare it with any target.

## 8. Recorded look-elsewhere data

Before target comparison, write:

- all 37 frozen curve nodes `(r,b,P,D,epsilon_D)`;
- the full momentum multiset and its sampled min/max;
- the number of resolved main-grid fold brackets;
- all refined bracket endpoints, widths and midpoint momentum estimates;
- even/odd discrepancies;
- branch margins and lapse ratios.

## 9. Mechanical outcomes

Assign the first applicable outcome:

1. `MOMENTUM_ENVELOPE_CONTROL_FAILED`;
2. `MOMENTUM_ENVELOPE_BASE_BRACKET_FAILED`;
3. `MOMENTUM_ENVELOPE_CURVE_NEWTON_OPEN`;
4. `MOMENTUM_ENVELOPE_BRANCH_OR_FULL_FAILED`;
5. `MOMENTUM_ENVELOPE_DERIVATIVE_OPEN`;
6. `MOMENTUM_ENVELOPE_FOLD_ENUMERATION_OPEN`;
7. `MOMENTUM_ENVELOPE_PARITY_FAILED`;
8. `MOMENTUM_ENVELOPE_ENUMERATED` only if every frozen node and fold gate
   passes.

## 10. Interpretation boundary

Outcome 8 is **DERIVED COMPUTATIONAL ON THE FROZEN DOMAIN**.  A resolved fold
is a genuine local obstruction/turning point of the canonical momentum map.

The finite grid plus sentinels is not an analytic global envelope on
`r in (-infinity,r1]`.  Absence of additional folds between samples or beyond
the sentinels remains **OPEN**; near-null convergence is only **PATTERN**.

Only the new targeted verifier is run.  The full suite is excluded.
