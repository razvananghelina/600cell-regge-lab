# Protocol: constrained-response auxiliary failure diagnostic

Date: 2026-08-21

Status: frozen before inspecting any finer stationarity derivative or
off-shell curvature ladder.

The first complete direct-action adversarial run is frozen in commit
`136de6d` with artifact

```text
reproducible/gravity_600cell_refined_h4_constrained_response_adversarial.json
  a23ef4cc23d08ad8768f1df66789aa900cdb95a7f3529486df80697a53b1fe81
docs/gravity/gravity_600cell_refined_h4_constrained_response_adversarial_first_result.md
  c4203c07b859ed323ee5049875d54d4894a1b815c3792b7c1d1de0e71677ad64
```

It passed the decisive matrix census but formally failed two auxiliary gates.
This diagnostic may explain or confirm those failures; it may not alter the
frozen `15/17` outcome.

## 1. Inputs and scope

Use exactly the scalar-action and geometry inputs frozen by the adversarial
protocol, plus the failed verifier itself at

```text
reproducible/verify_gravity_600cell_refined_h4_constrained_response_adversarial.py
  78f6b52f6f019a150a86ddadcb819b67c3757244c015687ab67f4649784ac53d
```

Reconstruct all 24 schedules and the same algebraic `20`-direction slice
`W=diag(P,Q)`.  Do not load the primary response matrices and do not construct
or compare any new response class.  The diagnostic concerns only (a) the ten
internal first derivatives at the static background and (b) the analytic
branch behaviour of the scalar action along the 210 directions used by the
failed audit.

## 2. Stationarity ladder

At both 180 and 220 decimal digits, for all ten raw internal log-coordinate
directions on all 24 schedules, compute centred first differences at

```text
h=(4e-15,2e-15,1e-15,5e-16,2.5e-16,1.25e-16).
```

Starting from `D_j`, apply

```text
R_j=(4 D_(j+1)-D_j)/3,
U_j=(16 R_(j+1)-R_j)/15,
V_j=(64 U_(j+1)-U_j)/63,
W_j=(256 V_(j+1)-V_j)/255.
```

Use `W_1` as the final derivative and freeze, componentwise,

```text
e_d=100 |W_1-W_0| + 1e-60 max(1,max_j |D_j|).
```

The stationary-zero hypothesis passes only if all 240 values satisfy

```text
|W_1| <= e_d
```

at both precisions and the two `W_1` values agree within the sum of their
envelopes.  This criterion does not merely ask for a decreasing sequence: a
resolved nonzero derivative converges away from zero and fails.

The `1e-60` floor is the already accepted upstream stationarity scale and is
also below the roughly 75 decimal digits stored for the selected masses.  A
smaller floor would claim information absent from the frozen inputs; it is
not evidence of exact symbolic stationarity.

Known controls:

- the same ladder applied to an even polynomial action must recover zero;
- adding a linear term `1e-20*x` must recover a resolved nonzero derivative
  and fail the stationary-zero gate.

## 3. Off-shell Lorentzian curvature diagnostic

At 180 decimal digits, retain every individual triangle curvature while
evaluating the complete scalar action.  For all 24 schedules and all 210
directions `e_a` and `e_a+e_b`, evaluate both signs at

```text
h=(1e-10,5e-11,2.5e-11).
```

Let `M_j` be the global maximum absolute imaginary individual curvature at
step `h_j`.  Let `M_base` be its maximum at the static background.  Require:

```text
M_base < 1e-100,
M_2 > 1e-12,
1.99 <= M_0/M_1 <= 2.01,
1.99 <= M_1/M_2 <= 2.01.                         (1)
```

For each schedule, direction and step, align triangle keys between the plus,
minus and base evaluations.  Define

```text
odd_scale=max_k |Im K_k(+h)-Im K_k(-h)|,
even_error=max_k |Im K_k(+h)+Im K_k(-h)-2 Im K_k(0)|.
```

Whenever `odd_scale>1e-100`, require at the finest step

```text
even_error/odd_scale < 1e-6.                     (2)
```

Across all evaluations also require relative imaginary complete action and
angle-identity residual below `1e-100`, and minimum angle argument above
`1e-2`.  Equations (1)--(2), together with a real complete action, identify a
smooth off-shell complex decomposition whose imaginary individual terms
vanish at the real background; they do not assert that every Lorentzian
deficit angle is itself an observable real number off shell.

Known controls:

- `K(x)=K0+i(7x+11x^3)` must pass the halving and oddness tests;
- `K_bad(x)=K0+i sign(x)` must fail the halving test by a resolved amount.

## 4. Frozen interpretation and outcomes

If both physical diagnostics and all known controls pass, the two failures of
the first audit are classified as:

- stationarity: finite-difference truncation from an under-resolved
  fourth-order derivative;
- curvature reality: an overstrong category error applied to individual
  off-shell terms even though the complete action stays real on a continuous
  safe branch.

Only then use

```text
REFINED_H4_CONSTRAINED_RESPONSE_AUXILIARY_FAILURES_RESOLVED.
```

If the known controls fail, use

```text
REFINED_H4_CONSTRAINED_RESPONSE_AUXILIARY_DIAGNOSTIC_CONTROL_FAILED.
```

Otherwise use

```text
REFINED_H4_CONSTRAINED_RESPONSE_AUXILIARY_FAILURES_CONFIRMED.
```

The last outcome leaves the direct adversarial result formally failed.  The
first outcome permits, but does not itself constitute, a separately
preregistered corrected adversarial run.  No full suite, root search, deferred
nonlinear census, nonhomogeneous spectrum or physical-constant extraction is
part of this diagnostic.
