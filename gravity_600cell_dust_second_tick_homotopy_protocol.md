# Preregistration: fixed canonical-target homotopy to the second dust tick

Date: 2026-08-16

Prior-art gate: `e760462`.

Accepted start: `46a7361`.

Direct-solver boundary: `6346ad0`.

Status: **frozen before evaluating `lambda>0`**.

## 1. Frozen homotopy

For each staircase parity load from the accepted first-tick artifact:

```text
a1 = log(L1/L0),
r1 = log(rho1/rho0),
t0 = the 30-component target momentum matched by the first tick,
t1 = mapped p_post,1, the desired second-tick incoming momentum.
```

Use exactly the 33 nodes

```text
lambda_j = j/32,  j=0,...,32.
```

At a fixed node define

```text
ell_j    = lambda_j*a1,
target_j = (1-lambda_j)*t0 + lambda_j*t1.
```

The two unknown absolute logarithms are

```text
b = log(L_upper/L0),
r = log(rho/rho0),
```

and the exact slab geometry is

```text
q_old       = exp(2*ell_j)*L0^2,
q_new       = exp(2*b)*L0^2,
pole        = exp(r)*rho0,
diagonal    = exp(ell_j+b)*L0^2-exp(r)*rho0.
```

The fixed equations are

```text
F0_j(b,r) = mean(g_pole[5]),
F1_j(b,r) = mean(p_pre[30]-target_j[30]).
```

The mass, action, maps and target endpoints are fixed.  Intermediate targets
are the componentwise linear interpolation above; no alternative path is
allowed.

## 2. Start and predictor

At `lambda_0=0` require the stored state `(b,r)=(a1,r1)` to reproduce all 35
internal equations and all 30 momentum equations within the endpoint gates.

For every `j>0`, use the accepted state at `j-1` unchanged as the sole seed.
Do not use the failed direct-solver trajectory, extrapolation, tangent
prediction, alternate seed, grid search, adaptive lambda, bisection or branch
restart.

## 3. Jacobian calibration

At every Newton state construct four central-difference matrices in `(b,r)`:

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Set

```text
epsilon = ||Jop1-Jop2||_2 + ||Jval1-Jval2||_2
          + ||Jop1-Jval1||_2 + 1e-60.
```

Require the existing entrywise factor-10 consistency test and
`s_min(Jop1)>100*epsilon`.  Every finite-difference evaluation must remain on
the complete Lorentzian/complex-angle branch.  Record all matrices, singular
values, determinant, condition number and error diagnostics.

## 4. Node corrector

At each `j>0`, apply the frozen Newton/Armijo rule

```text
Jop1*delta=-F,
alpha in (1,1/2,...,2^-10),
||F_trial||_infinity
  <= (1-alpha/4)*||F_current||_infinity.
```

Accept the first branch-valid damping endpoint.  Recalibrate after every
accepted correction.  Each node gets at most six accepted Newton iterations
and must reach `||F||_infinity<1e-25`.

No adaptive step size, coordinate rescaling, optimizer, Broyden update or
tolerance change is permitted.

## 5. Complete gate at every node

At every accepted node, not only at `lambda=1`, require

```text
max abs(30 diagonal equations) < 1e-60,
max abs(5 pole equations)      < 1e-25,
max abs(all 35 equations)      < 1e-25,
within-type residual spreads   < 1e-60,
||p_pre-target_j||_2           <= inherited first-tick junction bound,
spread(p_pre-target_j)         <= inherited first-tick junction bound,
calibrated endpoint rank       = 2.
```

Record at each node

```text
lambda, b, r,
L_upper/L_lower,
rho/rho1, tau/tau1,
all 35 residuals and 30 momentum residuals,
the endpoint Jacobian diagnostics.
```

## 6. Independent schedule gate

Run the complete 32-step homotopy separately for even and odd schedules.  At
every corresponding accepted node require

```text
abs(b_even-b_odd) < 1e-25,
abs(r_even-r_odd) < 1e-25.
```

At `lambda=1` additionally require pre/post momentum infinity differences
below `1e-22`.

## 7. Mechanical outcome hierarchy

Assign the first applicable outcome:

1. `SECOND_TICK_HOMOTOPY_CONTROL_FAILED` for a provenance, start, map or
   branch-control failure;
2. `SECOND_TICK_HOMOTOPY_JACOBIAN_OPEN` if any node Jacobian is unresolved;
3. `SECOND_TICK_HOMOTOPY_NEWTON_OPEN` if any fixed node misses its corrector
   gate;
4. `SECOND_TICK_HOMOTOPY_FULL_SUBSTITUTION_FAILED` if any reduced root fails
   a complete node gate;
5. `SECOND_TICK_HOMOTOPY_SCHEDULE_DEPENDENT` if the independent paths differ;
6. `SECOND_TICK_HOMOTOPY_STATIONARY` if all gates pass and endpoint
   `abs(b-a1)<=1e-20`;
7. `SECOND_TICK_HOMOTOPY_CONTINUED_CONTRACTION` if all gates pass and endpoint
   `b-a1<-1e-20`;
8. `SECOND_TICK_HOMOTOPY_TURNED_TO_EXPANSION` if all gates pass and endpoint
   `b-a1>1e-20`.

## 8. Interpretation boundary

A passing non-static endpoint is **DERIVED COMPUTATIONAL LOCAL** evidence for
a connected second canonical homothetic slab.  The intermediate homotopy is a
numerical device and must not be interpreted as 32 physical sub-ticks.

The result cannot establish global uniqueness, anisotropic stability,
refinement convergence or absolute/emergent time.  The initial `tau0` remains
external, and the selected relative lapse remains **STRUCTURAL / candidate
pseudo-constraint**.

Only the new targeted verifier is run.  The full suite remains excluded.
