# Preregistration: weak-lapse jet of the three-step canonical map

Date: 2026-08-16

Prior-art gate: `90757a6`.

Status: frozen before evaluating any scaled-lapse non-static slab.

## 1. Frozen inputs and control family

Require exact SHA-256 values

```text
regular-lapse identity:
5079428fade247f730ebc07e5e2eae388b48045cd5201e84afb3186bfc248a51,

accepted tick 1:
4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9,

accepted tick 2:
936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70,

accepted tick 3:
ebf2f1a11b9a4e9c76fb1ce33066c0782429cf6500770df7bbe4d92de4a050c0,

gluing map:
a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77.
```

Use exactly

```text
lambda in {1/2,1/4,1/8}.
```

For each `lambda`, set

```text
rho_base=lambda^2*rho0,
tau_base=lambda*tau0,
k_lambda=lambda*epsilon3*L0*tau0/4.
```

Keep `L0`, the dust mass and every carrier/action convention fixed.

## 2. Exact scaled static controls

For both parities and every `lambda`, freshly evaluate

```text
lower log = upper log = 0,
relative lapse log = 0.
```

Require all 35 internal residuals below `1e-60`, Lorentzian branch gates and

```text
pre momentum  = -k_lambda,
post momentum = +k_lambda
```

componentwise to relative error below `1e-60`.

## 3. Sequential canonical map

At step `n=1,2,3`, given lower log `a_(n-1)` and incoming 30-vector `p_n`,
solve for absolute upper log `a_n` and relative lapse log `r_n` using

```text
q_old    = exp(2*a_(n-1))*L0^2,
q_new    = exp(2*a_n)*L0^2,
rho_n    = lambda^2*rho0*exp(r_n),
diagonal = exp(a_(n-1)+a_n)*L0^2-rho_n,

G(a_n,r_n)=mean(five pole equations)=0,
P_pre(a_n,r_n)=mean(p_n).
```

The input for step 1 is the exact static output `+k_lambda`.  Inputs for steps
2 and 3 are the preceding freshly computed post-momenta after applying the
independent old-to-final orbit map.  Thus no momentum target is fitted or
imported from the `lambda=1` trajectory.

## 4. Frozen continuation seeds and solver

Use the accepted `lambda=1` absolute states only as zero-order branch seeds:

```text
a_seed(n,lambda) = lambda^2*a_n(lambda=1),
r_seed(n,lambda) = lambda^2*r_n(lambda=1).
```

No alternate seed, stationary-root search, restart or solver endpoint is
allowed.  This is a connected-branch asymptotic control, not blind branch
discovery.

Use the same four derivative steps as the accepted ticks:

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Require the standard entrywise derivative agreement and
`s_min>100*epsilon`.  Permit at most eight accepted Newton iterations per
step.  Try Armijo factors `1,1/2,...,2^-10` in order, using the frozen decrease
rule

```text
||F_trial||_infinity <= (1-alpha/4)*||F_current||_infinity.
```

Stop only at `||F||_infinity<1e-25`.

## 5. Full gates

For all 18 solves (three lambdas, three steps, two parities) require

```text
max abs(30 diagonal equations) < 1e-60,
max abs(5 pole equations)      < 1e-25,
max abs(all 35 equations)      < 1e-25,
diagonal spread                < 1e-60,
pole spread                    < 1e-60,
norm_2(pre-target)             < 1e-24,
spread(pre-target)             < 1e-24.
```

Every action, derivative and trial state must retain the Lorentzian branch.
For each `(lambda,n)`, require even/odd agreement

```text
|a_even-a_odd| < 1e-25,
|r_even-r_odd| < 1e-25,
max abs(pre_even-pre_odd)  < 1e-22,
max abs(post_even-post_odd)< 1e-22.
```

## 6. Frozen recurrence observables

For each parity and `lambda`, set `a_0=r_0=0` and

```text
u_n=a_n-a_(n-1),
v_n=r_n-r_(n-1).
```

Record the following eleven nontrivial ratios and targets:

```text
u_2/u_1 -> 2,       u_3/u_1 -> 3,
a_2/u_1 -> 3,       a_3/u_1 -> 6,
v_2/v_1 -> 3,       v_3/v_1 -> 5,
r_2/v_1 -> 4,       r_3/v_1 -> 9,
p_post,1/k_lambda -> 3,
p_post,2/k_lambda -> 5,
p_post,3/k_lambda -> 7.
```

Also record the scaled leading coefficients `u_1/lambda^2` and
`v_1/lambda^2` without assigning them a target.

## 7. Quadratic asymptotic test

For each ratio `Q(lambda)` with integer target `q`, define

```text
e(lambda)=abs(Q(lambda)-q),
order_12=log2(e(1/2)/e(1/4)),
order_23=log2(e(1/4)/e(1/8)).
```

Errors are resolved if they exceed `1e-40`.  A ratio has a quadratic integer
limit only if both errors are resolved, strictly decrease under both halvings
and

```text
1.8 <= order_12 <= 2.2,
1.8 <= order_23 <= 2.2.
```

Compute two target-independent Richardson intercepts

```text
Q0_coarse=(4*Q(1/4)-Q(1/2))/3,
Q0_fine  =(4*Q(1/8)-Q(1/4))/3,
epsilon_R=abs(Q0_fine-Q0_coarse)+1e-40.
```

The frozen integer is consistent with the extrapolated intercept only if

```text
abs(Q0_fine-q) <= 10*epsilon_R.
```

This is an empirical asymptotic certificate, not a rigorous error bound.

## 8. Mechanical outcomes

Assign exactly one:

1. `WEAK_LAPSE_RECURRENCE_CONTROL_FAILED`;
2. `WEAK_LAPSE_RECURRENCE_JACOBIAN_OPEN`;
3. `WEAK_LAPSE_RECURRENCE_NEWTON_OPEN`;
4. `WEAK_LAPSE_RECURRENCE_FULL_GATE_FAILED`;
5. `WEAK_LAPSE_RECURRENCE_SCHEDULE_DEPENDENT`;
6. `WEAK_LAPSE_INTEGER_TREND_ONLY` if all final errors are smaller than the
   initial errors but any quadratic/order/Richardson gate fails;
7. `WEAK_LAPSE_QUADRATIC_INTEGER_LAW` if all eleven ratios pass every frozen
   asymptotic gate.

Outcome 7 establishes, computationally on this branch,

```text
p_post,n/k_lambda = 2n+1+O(lambda^2),
u_n/u_1           = n+O(lambda^2),
v_n/v_1           = 2n-1+O(lambda^2),
```

and hence triangular/square cumulative logs to leading weak-lapse order.  It
does not prove the exact all-order recurrence, spatial refinement, Einstein
continuum convergence or emergent time.

Only the new targeted verifier will be run.  The full suite will not be run.
