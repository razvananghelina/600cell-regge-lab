# Preregistration: local canonical correction of the third dust tick

Date: 2026-08-16

Prior-art gate: `7b9a676`.

Target-independent roots: `3401137`.

Disclosed target comparison: `7cf4e27`.

Status: **cleanly frozen before evaluating any target-corrected third-slab
state**.

## Exact carrier, conserved mass and state

Use both derived order-24 staircase schedules and the complete 100-decimal
Lorentzian Regge plus dust action.  Hold fixed the lower boundary

```text
B2 = log(L2/L0)
   = -9.34818705890582713633822299265753373027428194008991504419612e-6
```

and the original conserved mass

```text
M=(90/pi)*epsilon3*L0.
```

The sole unknowns are absolute logarithms

```text
C = log(L3/L0),
R = log(rho3/rho0).
```

Use exactly

```text
q_old    = exp(2*B2)*L0^2,
q_new    = exp(2*C)*L0^2,
rho3     = exp(R)*rho0,
diagonal = exp(B2+C)*L0^2-rho3.
```

## Sole seed and target

Require input hashes

```text
stationary roots:
02d4589a7df0851c67a31fc0a41c5ef8851a82c758214c1c5e8729afddfe479f,

accepted second tick:
936984bc84a714140ce16917ee559b346b3c0d4a5ba92d8fb723398a120f8e70,

gluing map:
a5a22d219b71e49c154c1ef80ed9da93b1aef0b93cd2d6ed22f041b71f62db77,

target comparison:
4d1f81dafcab9d3aa40ff08fdaaad90b80235809dd32becd790bdee1704ab6cf.
```

Use only committed contracting root index 0:

```text
C_seed = -1.86964989240472221703406162566722667346751807829667980709325e-5,
R_seed = -1.42370275520098029961300545242474815338378370661665379256974e-5.
```

It is selected by `C_seed<B2`, not by target proximity.  No time-reversal
root, failed endpoint, alternative bracket or restart is permitted.

For each parity construct the 30-component target by applying the independent
old-to-final map to the accepted second tick's post-momentum.

## Reduced equations and frozen Newton

Define

```text
G(C,R) = mean of five pole equations,
P(C,R) = mean of thirty pre-momenta,
F(C,R) = (G(C,R), P(C,R)-mean(target)).
```

Drive Newton only with `F`, but accept only after substituting all 35 internal
equations and all 30 momentum components.

Use the four central-difference steps

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Apply the identical entrywise consistency and singular-value error gates as
the accepted first two ticks.  Require the smaller singular value to exceed
`100*epsilon` at every attempted and endpoint Jacobian.

Allow at most eight accepted iterations.  Try Armijo dampings
`1,1/2,...,2^-10` in order and accept the first with

```text
||F_trial||_infinity <= (1-alpha/4)*||F_current||_infinity.
```

Stop successfully only at `||F||_infinity<1e-25`.  No parameter may be
retuned after evaluation.

## Full gates

At each parity endpoint require

```text
max abs(30 diagonal equations) < 1e-60,
max abs(5 pole equations)      < 1e-25,
max abs(all 35 equations)      < 1e-25,
diagonal spread                < 1e-60,
pole spread                    < 1e-60,
norm_2(pre-target)             <= inherited junction_bound,
spread(pre-target)             <= inherited junction_bound.
```

Every action/derivative/trial state must retain all Lorentzian branch gates.
Require parity agreement

```text
|C_even-C_odd| < 1e-25,
|R_even-R_odd| < 1e-25,
max abs(pre_even-pre_odd)   < 1e-22,
max abs(post_even-post_odd) < 1e-22.
```

Report

```text
u3=C-B2, v3=R-R2,
L3/L2=exp(u3), tau3/tau2=exp(v3/2),
u3/A1, v3/R1, C/A1, R/R1.
```

## Predictions frozen before new evaluation

The linear correction computed solely from committed artifacts is

```text
delta_C approximately +9.981943499e-11,
delta_R approximately -1.779618554e-5.
```

The two-tick pattern predicts diagnostically

```text
u3/A1 approximately 3,
v3/R1 approximately 5,
C/A1 approximately 6,
R/R1 approximately 9.
```

No agreement is required for acceptance.  Failure kills the integer-sequence
pattern but does not override the canonical gates.

## Mechanical outcomes

Assign exactly one:

1. `THIRD_TICK_LOCAL_CONTROL_FAILED`;
2. `THIRD_TICK_LOCAL_JACOBIAN_OPEN`;
3. `THIRD_TICK_LOCAL_NEWTON_OPEN`;
4. `THIRD_TICK_LOCAL_FULL_SUBSTITUTION_FAILED`;
5. `THIRD_TICK_LOCAL_SCHEDULE_DEPENDENT`;
6. `THIRD_HOMOTHETIC_TICK_ACCEPTED`.

Outcome 6 establishes a third consecutive homogeneous canonical slab only.
General recurrence, continuum/refinement stability, anisotropic modes and an
absolute clock remain **OPEN**.

Only the new targeted verifier will be run.  The full suite will not be run.
