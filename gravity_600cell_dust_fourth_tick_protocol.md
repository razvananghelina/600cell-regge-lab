# Preregistration: out-of-sample fourth weak-lapse dust tick

Date: 2026-08-16

Prior-art gate: `40d77c7`.

Status: frozen before evaluating any fourth non-static slab.

## 1. Frozen inputs

Require the committed weak-lapse artifact

```text
reproducible/gravity_600cell_dust_weak_lapse_recurrence.json
SHA-256 = 500be1c4e2d7ec4104b9773bc1cfc71065c9d930607119eb616367d18fa5d8f9,
outcome  = WEAK_LAPSE_QUADRATIC_INTEGER_LAW,
tests    = passed = 5,
tick4_target_parsed = false.
```

Freeze the same complete-action evaluator, old-to-final orbit maps, precision,
derivative steps and lambdas used there:

```text
base solver SHA-256
  cef59fa0bc3a1c8fa3be0193234371b7dda303a0ec72683ddcdd88bcb40f3725,

gluing artifact SHA-256
  a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77,

lambda in {1/2,1/4,1/8},
DPS = 100.
```

Keep `L0`, the dust mass, carrier, action and branch conventions fixed.

## 2. Fourth canonical seam

For every `lambda` and parity, parse the three committed states.  Re-evaluate
the third slab from its committed second and third states and require its
pre/post momenta to reproduce the stored 50-digit vectors to infinity error
below `1e-45`.

Set

```text
lower log a = committed a3,
incoming target = committed p_post,3 after the independent
                  old-to-final orbit map.
```

Solve only for the fourth absolute log scale `a4` and fourth absolute
relative-lapse log `r4`:

```text
q_old    = exp(2*a3)*L0^2,
q_new    = exp(2*a4)*L0^2,
rho4     = lambda^2*rho0*exp(r4),
diagonal = exp(a3+a4)*L0^2-rho4,

G(a4,r4)     = mean(five pole equations) = 0,
P_pre(a4,r4) = mean(incoming target).
```

No fixed-lapse root and no `lambda=1` fourth state may be evaluated.

## 3. Frozen seed and solver

For each parity and lambda, obtain `u1=a1` and `v1=r1` from the committed
first scaled solve.  Use the leading recurrence only as the unique branch
seed:

```text
a4_seed = 10*u1,
r4_seed = 16*v1.
```

No alternative seed, root enumeration, restart or endpoint selection is
allowed.

Use exactly

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Require the standard entrywise derivative agreement and
`s_min>100*epsilon`.  Permit at most eight accepted Newton iterations.  Try
Armijo factors `1,1/2,...,2^-10` in order with

```text
||F_trial||_infinity <= (1-alpha/4)*||F_current||_infinity.
```

Stop only at `||F||_infinity<1e-25`.

## 4. Full gates

For all six new solves require

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
At each lambda require even/odd agreement

```text
|a4_even-a4_odd|              < 1e-25,
|r4_even-r4_odd|              < 1e-25,
max abs(pre4_even-pre4_odd)   < 1e-22,
max abs(post4_even-post4_odd) < 1e-22.
```

## 5. Frozen out-of-sample predictions

For each parity and lambda form

```text
u4 = a4-a3,
v4 = r4-r3,
k  = lambda*epsilon3*L0*tau0/4.
```

The five predictions, derived only from the committed `n<=3` law, are

```text
u4/u1       -> 4,
a4/u1       -> 10,
v4/v1       -> 7,
r4/v1       -> 16,
p_post,4/k  -> 9.
```

## 6. Frozen quadratic and uncertainty rules

For each ratio `Q(lambda)` and frozen integer `q`, define

```text
e(lambda)=abs(Q(lambda)-q),
order_12=log2(e(1/2)/e(1/4)),
order_23=log2(e(1/4)/e(1/8)).
```

Require all errors above `1e-40`, strict decrease under both halvings, and

```text
1.8 <= order_12 <= 2.2,
1.8 <= order_23 <= 2.2.
```

Compute

```text
Q0_coarse=(4*Q(1/4)-Q(1/2))/3,
Q0_fine  =(4*Q(1/8)-Q(1/4))/3,
epsilon4 =abs(Q0_fine-Q0_coarse)+1e-40.
```

Require both an internal convergence check

```text
abs(Q0_fine-q) <= 10*epsilon4
```

and the fixed external prediction band

```text
B_train = 4.6222921056804246599831556548181231e-10,
abs(Q0_fine-q) <= B_train.
```

`B_train` is ten times the largest Richardson remainder among all eleven
committed `n<=3` observables.  It is fixed here before tick four and must not
be enlarged after evaluation.  These are empirical asymptotic controls, not
rigorous truncation-error bounds.

## 7. Mechanical outcomes

Assign exactly one:

1. `FOURTH_TICK_CONTROL_FAILED`;
2. `FOURTH_TICK_JACOBIAN_OPEN`;
3. `FOURTH_TICK_NEWTON_OPEN`;
4. `FOURTH_TICK_FULL_GATE_FAILED`;
5. `FOURTH_TICK_SCHEDULE_DEPENDENT`;
6. `FOURTH_TICK_INTEGER_TREND_ONLY` if every error decreases but any frozen
   order or Richardson gate fails;
7. `FOURTH_TICK_WEAK_LAPSE_PREDICTION_REFUTED` otherwise;
8. `FOURTH_TICK_WEAK_LAPSE_PREDICTION_CONFIRMED` only if all five ratios pass
   every frozen trend, order, internal Richardson and external-band gate.

Outcome 8 validates the leading recurrence at one iteration index not used to
identify it:

```text
p_post,4/k = 9+O(lambda^2),
u4/u1      = 4+O(lambda^2),
v4/v1      = 7+O(lambda^2).
```

It does not prove an exact or infinite recurrence, spatial refinement,
continuum Einstein dynamics or emergent time.

Only the new targeted verifier will be run.  The full suite will not be run.
