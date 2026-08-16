# Preregistration: second canonical homothetic 600-cell dust tick

Date: 2026-08-16

Prior-art gate: `1865e13`.

Upstream first-tick result: `46a7361`.

Status: **frozen before evaluating the second-tick seed or equations**.

## 1. Fixed input, variables and geometry

For each derived staircase parity load the first-tick values `(a1,r1)` and its
30-component post-momentum from the committed accepted artifact.  Load the
vertex-derived old-to-final orbit map from the committed two-slab control and
define

```text
target[old] = p_post,1[old_to_final[old]].
```

The sole unknown pair is

```text
u = log(L2/L1),
v = log(rho2/rho1).
```

Use exactly

```text
q_old       = exp(2*a1)*L0^2,
q_new       = exp(2*(a1+u))*L0^2,
pole        = exp(r1+v)*rho0,
diagonal    = exp(2*a1+u)*L0^2-exp(r1+v)*rho0.
```

The fixed reduced equations are

```text
F0(u,v) = mean(g_pole[5]),
F1(u,v) = mean(p_pre,2[30]-target[30]).
```

No mass, action coefficient, lower length, incoming momentum, orbit map,
internal length, lapse prior or target coefficient may vary.

## 2. Sole seed

Use only

```text
u_seed = a1,
v_seed = 0.
```

This repeats the preceding logarithmic scale increment while carrying forward
the preceding proper tick duration.  It is a deterministic constant-velocity,
constant-step predictor, not a fitted estimate.

No alternate seed, grid search, bracket, random restart or branch search is
allowed.  Failure from this seed is reported as open for this local
continuation; it is not a global no-go theorem.

## 3. Upstream and reconstruction controls

Before evaluating the second-tick seed require:

1. the first-tick artifact has outcome
   `HOMOTHETIC_CANONICAL_LAPSE_SELECTED`, `7/7`, and the expected provenance;
2. the two-slab gluing artifact passes and supplies a bijective old-to-final
   map for both parities;
3. a generalized evaluator at lower scale `L0`, upper scale `L1` and pole
   `rho1` reproduces the committed first tick with
   `max|35 residuals|<1e-24` and each pre/post momentum within `1e-45` of its
   stored 50-digit value;
4. reconstructed even/odd first-tick data agree within `1e-24`;
5. the 30-component second target is finite and complete.

Any failure assigns the control outcome and forbids a physical conclusion.

## 4. Frozen Jacobian calibration

At every Newton state build central-difference matrices in `(u,v)` with the
already validated dimensionless steps

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Let

```text
d_op    = J_op_primary-J_op_shadow,
d_val   = J_val_primary-J_val_shadow,
d_cross = J_op_primary-J_val_primary,
epsilon = ||d_op||_2+||d_val||_2+||d_cross||_2+1e-60.
```

Require each entry of `d_cross` to be no larger than ten times the sum of the
corresponding absolute `d_op` and `d_val` entries plus `1e-60`.  Invert only
when the smallest singular value of `J_op_primary` is greater than
`100*epsilon`.

Every perturbed evaluation must pass the full Lorentzian/complex-angle branch
gate.  Serialize all four matrices, singular values, determinant, condition
number, branch diagnostics and `epsilon` at every attempt.

## 5. Frozen Newton solve

At each accepted state solve

```text
J_op_primary*delta = -F.
```

Try damping factors `1,1/2,...,2^-10` in that order and accept the first
branch-valid endpoint satisfying

```text
||F_trial||_infinity
  <= (1-alpha/4)*||F_current||_infinity.
```

Recalibrate the complete four-matrix Jacobian after every accepted step.
Stop successfully at `||F||_infinity<1e-25`.  Otherwise stop after eight
accepted iterations, complete damping failure, a branch failure or an
unresolved Jacobian.

No coordinate rescaling, optimizer, Broyden update, tolerance change or
manual correction is permitted.

## 6. Complete endpoint substitution

At a reduced root require

```text
max abs(30 diagonal equations) < 1e-60,
max abs(5 pole equations)      < 1e-25,
max abs(all 35 equations)      < 1e-25,
within-type residual spreads   < 1e-60,
||p_pre,2-target||_2           <= inherited first-tick junction bound,
spread(p_pre,2-target)         <= inherited first-tick junction bound.
```

The endpoint Jacobian must again have calibrated rank two.  Record

```text
u, L2/L1, L2/L0,
v, rho2/rho1, rho2/rho0, tau2/tau1, tau2/tau0,
all 35 internal residuals,
all 30 pre/post momenta,
all 30 junction residuals.
```

## 7. Parity gate

Solve even and odd schedules separately.  Require

```text
abs(u_even-u_odd)             < 1e-25,
abs(v_even-v_odd)             < 1e-25,
max abs(p_pre_even-p_pre_odd) < 1e-22,
max abs(p_post_even-p_post_odd)<1e-22.
```

## 8. Mechanical outcome hierarchy

Assign exactly one outcome in this order:

1. `SECOND_TICK_CONTROL_FAILED` for any provenance, reconstruction, target or
   branch-control failure;
2. `SECOND_TICK_JACOBIAN_OPEN` if any calibrated Jacobian is not rank two;
3. `SECOND_TICK_NEWTON_OPEN` if the deterministic solve misses its residual
   gate;
4. `SECOND_TICK_FULL_SUBSTITUTION_FAILED` if a reduced root fails a complete
   65-component gate;
5. `SECOND_TICK_SCHEDULE_DEPENDENT` if the two endpoints fail the parity gate;
6. `SECOND_TICK_STATIONARY` if every gate passes and `abs(u)<=1e-20`;
7. `SECOND_TICK_CONTINUED_CONTRACTION` if every gate passes and `u<-1e-20`;
8. `SECOND_TICK_TURNED_TO_EXPANSION` if every gate passes and `u>1e-20`.

The sign thresholds are display/classification gates, not solve criteria.

## 9. Interpretation boundary

Outcomes 7 or 8 establish a **DERIVED COMPUTATIONAL LOCAL** second non-static
canonical move in the homothetic subspace.  Outcome 6 establishes an accepted
stationary second move, not continuing dynamics.  Outcomes 1--5 leave local
iteration open or fail it under the frozen predictor.

Even a passing non-static outcome does not establish global uniqueness,
anisotropic stability, refinement convergence, an absolute time unit or
emergent time.  The selected lapse remains **STRUCTURAL / candidate
pseudo-constraint** until its weak singular scale is tracked under further
iteration and refinement.

Only the new registered verifier will be run.  The full suite is excluded by
the user's instruction.
