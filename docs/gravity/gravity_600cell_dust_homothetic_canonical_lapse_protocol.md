# Preregistration: two-variable homothetic canonical-lapse solve

Date: 2026-08-16

Prior-art gate: `c7f3e29`.

Upstream fixed-lapse result: `b788258`.

Status: **frozen before any evaluation with `rho_next != rho0`**.

## 1. Fixed equation, seed and target

For each derived schedule parity retain exactly the complete 100-decimal
Lorentzian Regge+dust evaluator and homothetic geometry of the upstream
verifiers.  The unknowns are

```text
y = (s,z),
s = log(L_+/L0),
z = log(rho_next/rho0),
rho0=(0.0102)^2.
```

Set

```text
q_old       = L0^2,
q_new       = exp(2s)*L0^2,
pole        = rho0*exp(z),
diagonal    = exp(s)*L0^2-rho0*exp(z).
```

Load the 30-component forward target exactly as in the fixed-lapse root
verifier from the committed two-slab artifact.  Define

```text
F0(y) = mean(g_pole[5]),
F1(y) = mean(p_pre[30]-p_target[30]).
```

Use the sole seed

```text
s = the committed fixed-lapse root,
z = 0.
```

No alternate seed, bracket, mass, target coefficient, lapse prior or
internal length is allowed.

## 2. Base reproduction control

Before differentiating in `z`, require at the seed:

- maximum internal residual below `1e-25`;
- the fresh `s` agrees with the committed root string exactly as loaded;
- the fresh uniform momentum mismatch agrees with the committed
  per-component mismatch within `1e-20`;
- both schedules remain Lorentzian and reproduce each other within `1e-24`;
- all target maps and uncertainty norms are loaded from passing artifacts.

Stop with a control failure if any item fails.

## 3. Frozen Jacobian calibration

For a central-difference matrix in `(s,z)`, use these coordinate-wise steps:

```text
operational primary = 5e-9,
operational shadow  = 1e-8,
validation primary  = 1.5e-8,
validation shadow   = 3e-8.
```

At every Newton iterate evaluate all four matrices.  Let

```text
d_op    = J_op_primary-J_op_shadow,
d_val   = J_val_primary-J_val_shadow,
d_cross = J_op_primary-J_val_primary,
epsilon = ||d_op||_2+||d_val||_2+||d_cross||_2+1e-60.
```

Require every entry of `d_cross` to be at most ten times the sum of the
corresponding absolute `d_op`, `d_val` entries plus `1e-60`.  Invert only if
the smallest singular value of `J_op_primary` exceeds `100*epsilon`.

Every perturbed evaluation must pass the full Lorentzian and complex-angle
branch gates.  Report both singular values, determinant, condition number,
all four matrices and `epsilon` at every accepted iterate.

## 4. Frozen Newton corrector

At each iterate solve

```text
J_op_primary * delta = -F.
```

Try damping factors `1,1/2,...,2^-10` in that order.  Accept the first
branch-valid trial satisfying

```text
||F_trial||_infinity
  <= (1-alpha/4)*||F_current||_infinity.
```

Recalibrate the full four-matrix Jacobian after every accepted step.  Stop
successfully when `||F||_infinity<1e-25`, otherwise stop after eight accepted
Newton iterations, failure of all damping values, a branch failure or a
Jacobian error-band failure.

No SciPy optimizer, Broyden update, random restart, coordinate rescaling or
tolerance change is permitted.

## 5. Complete endpoint substitution

At an accepted endpoint require:

```text
max abs(30 diagonal equations) < 1e-60,
max abs(5 pole equations)      < 1e-25,
max abs(all 35 equations)      < 1e-25,
within-type residual spreads   < 1e-60,
||p_pre-p_target||_2
    <= 10*cusp_uncertainty_norm,
spread(p_pre-p_target)
    <= 10*cusp_uncertainty_norm.
```

The endpoint Jacobian must again have calibrated rank two.  Define the lapse
shift as resolved if `abs(z)>1e-20`; this fixed display gate is not used to
force convergence.

Record

```text
s,
L_+/L0=exp(s),
z,
rho_next/rho0=exp(z),
tau_next/tau0=exp(z/2),
all 35 internal residuals,
all 30 pre/post momenta,
all 30 junction residuals.
```

## 6. Parity gate

Require both schedules to converge and

```text
abs(s_even-s_odd)                 < 1e-25,
abs(z_even-z_odd)                 < 1e-25,
max abs(p_pre_even-p_pre_odd)     < 1e-22,
max abs(p_post_even-p_post_odd)   < 1e-22.
```

## 7. Mechanical outcome hierarchy

Assign exactly one outcome in this order:

1. `CANONICAL_LAPSE_CONTROL_FAILED` if an upstream, target, reproduction or
   branch control fails;
2. `CANONICAL_LAPSE_JACOBIAN_OPEN` if either calibrated Jacobian is not
   rank two;
3. `CANONICAL_LAPSE_NEWTON_OPEN` if the deterministic solve does not meet
   its reduced residual gate;
4. `CANONICAL_LAPSE_FULL_SUBSTITUTION_FAILED` if a reduced root fails any
   complete 65-component gate;
5. `CANONICAL_LAPSE_SCHEDULE_DEPENDENT` if the separate endpoints fail the
   parity gate;
6. `CANONICAL_ROOT_LAPSE_SHIFT_UNRESOLVED` if every physical equation passes
   but `abs(z)<=1e-20`;
7. `HOMOTHETIC_CANONICAL_LAPSE_SELECTED` only if every gate passes and the
   lapse shift is resolved.

## 8. Interpretation boundary

Outcome 7 establishes a **DERIVED COMPUTATIONAL LOCAL** rank-two canonical
root in the homothetic subspace.  It supplies the first canonically glued
non-static next slab and a next-lapse ratio selected by discrete consistency.

It does not establish full 65-variable global uniqueness, refinement
stability, an exact continuum gauge interpretation, a fundamental initial
tick or an absolute time unit.  The selected lapse must be labelled
**STRUCTURAL / candidate pseudo-constraint** until it survives those tests.

Only the new targeted verifier will be run.  It must be registered exactly
once; the full suite will not be run.
