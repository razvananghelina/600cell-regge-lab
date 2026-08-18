# Preregistration: local correction of the contracting second-tick root

Date: 2026-08-16

Prior-art and mass-framing gate: `fcc4d7c`.

Frozen stationary-target comparison: `1d32334`, outcome
`STATIONARY_SECOND_TICK_NO_HIT`.

Status: frozen before evaluating any target-corrected state.

## 1. Exact carrier, mass and unknowns

Use both derived order-24 staircase schedules, the complete 100-decimal
Lorentzian Regge plus dust action and the accepted first-tick geometry.  Hold
fixed:

```text
a1 = log(L1/L0),
M  = (90/pi)*epsilon3*L0.
```

The mass is conserved and is not recomputed from `L1` or `L2`.

The sole two unknowns are absolute logarithms

```text
b = log(L2/L0),
r = log(rho2/rho0).
```

For each state use exactly

```text
q_old    = exp(2*a1)*L0^2,
q_new    = exp(2*b)*L0^2,
rho2     = exp(r)*rho0,
diagonal = exp(a1+b)*L0^2-rho2.
```

No other edge, mass, target or coefficient is varied.

## 2. Sole seed and non-cherry-picking rule

Load the target-independent root artifact with SHA-256
`0ec5ba520ea25b39dd6cfd3c349d49fe480df2abee359854e1316b5af4d9fa2f`.
Use root index 0 and the inherited lapse as the sole seed:

```text
b_seed = -9.348231422816359925124508959895758699419300779844491278e-6,
r_seed = -3.55925313517063343725030533963917396571974345422547402551491e-6.
```

This root was selected before target comparison by `b_seed<a1`, hence by
contracting orientation.  Root 1 is the time-reversal branch and is not an
alternate seed.  No failed direct-solve or homotopy endpoint may be used.

## 3. Canonical target and reduced equations

For each parity map the accepted first tick's complete 30-component
`post_momentum` with the gluing artifact's `old_to_final_orbit_map`:

```text
target[i] = first_post[map[i]].
```

Define

```text
G(b,r) = mean of the five pole equations,
P(b,r) = mean of the thirty pre-momenta,
F(b,r) = (G(b,r), P(b,r)-mean(target)).
```

The two reduced equations only drive Newton.  Acceptance is always checked
on all 35 internal equations and all 30 momentum components.

## 4. Frozen derivative and Newton rules

For both coordinates use the four central-difference steps

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Use the identical entrywise consistency test and singular-value error envelope
of the accepted first-tick solver.  Require the smaller singular value to
exceed `100*epsilon` at every attempted and endpoint Jacobian.

At most eight accepted Newton iterations are allowed.  At each iteration try
Armijo dampings `1,1/2,...,2^-10` in that order and accept the first satisfying

```text
||F_trial||_infinity <= (1-alpha/4)*||F_current||_infinity.
```

Stop successfully at `||F||_infinity<1e-25`.  No restart, alternative seed,
extra iteration, derivative step or optimizer is allowed.

From only the already committed root, Jacobian and target mismatch, the
first linear correction predicts diagnostically

```text
delta_b approximately +4.4363871e-11,
delta_r approximately -1.0677721e-5.
```

This prediction is not an acceptance target.

## 5. Full acceptance gates

At the endpoint require separately, for each parity:

```text
max abs(30 diagonal equations) < 1e-60,
max abs(5 pole equations)      < 1e-25,
max abs(all 35 equations)      < 1e-25,
spread of diagonal equations   < 1e-60,
spread of pole equations       < 1e-60,
norm_2(pre-target)             <= inherited first-tick junction_bound,
spread(pre-target)             <= inherited first-tick junction_bound.
```

Every evaluated state and derivative displacement must retain the Lorentzian
branch and angle-cut gates.  The two parity solutions must agree within

```text
|b_even-b_odd| < 1e-25,
|r_even-r_odd| < 1e-25,
max abs(pre_even-pre_odd)  < 1e-22,
max abs(post_even-post_odd)< 1e-22.
```

Report the relative tick quantities

```text
u2=b-a1,  v2=r-r1,
L2/L1=exp(u2),  tau2/tau1=exp(v2/2),
u2/a1.
```

The last ratio is diagnostic, not fitted.

## 6. Mechanical outcome hierarchy

Assign exactly one, in order:

1. `SECOND_TICK_LOCAL_CONTROL_FAILED`;
2. `SECOND_TICK_LOCAL_JACOBIAN_OPEN`;
3. `SECOND_TICK_LOCAL_NEWTON_OPEN`;
4. `SECOND_TICK_LOCAL_FULL_SUBSTITUTION_FAILED`;
5. `SECOND_TICK_LOCAL_SCHEDULE_DEPENDENT`;
6. `SECOND_HOMOTHETIC_TICK_ACCEPTED`.

Outcome 6 establishes one locally selected second homogeneous tick at fixed
conserved mass.  It does not establish a global branch, arbitrary-tick
evolution, refinement stability, an absolute clock, `c` or Planck units.

Only the new targeted verifier will be run.  The full suite will not be run.
