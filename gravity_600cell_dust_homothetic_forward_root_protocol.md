# Preregistration: homothetic forward-root and canonical-junction test

Date: 2026-08-16

Prior-art gate: `2ce58bd`.

Upstream local-response result: `086009a`.

Status: **frozen before evaluating the newly predicted nonzero root
bracket**.

## 1. Fixed carrier, action and target

Use both derived order-24 staircase schedules and exactly the 100-decimal
Lorentzian Regge plus fixed-mass dust action, angle branch, orbit ordering
and logarithmic gradients used in the committed local-response verifier.

Hold fixed

```text
L_- = L0,
tau = 0.0102,
rho = tau^2,
M   = (90/pi)*(2*pi-5*acos(1/3))*L0.
```

For the sole unknown `s=log(L_+/L0)`, impose the exact geometric ansatz

```text
q_old       = L0^2,
q_new       = exp(2s)*L0^2,
pole        = rho,
diagonal    = exp(s)*L0^2-rho.
```

Let `E(s)` be the globally restricted lapse equation already defined by the
committed protocol.  The candidate root is required to satisfy the complete
35-vector `g_internal(s)`, not just `E(s)`.

Load `p_post(static)`, its uncertainty norm and the derived
`old_to_final_orbit_map` from
`reproducible/gravity_600cell_dust_two_slab_gluing.json`.  Construct the
forward target in the new old-boundary ordering as

```text
p_target[i] = p_post_static[old_to_final_orbit_map[i]].
```

No target component, scale coefficient, mass, lapse or internal edge is
optimized.

## 2. Frozen nonzero bracket and solver

For each parity use exactly

```text
s_left  = -1/40000  = -2.5e-5,
s_right = -1/640000 = -1.5625e-6.
```

The left point was part of the committed derivative protocol.  The right
point has not been evaluated.  This bracket was chosen from the disclosed
quadratic prediction `s approximately -3.116e-6` and is therefore labelled
**PATTERN-SEEDED**, not blind.

Evaluate both endpoints afresh.  Continue only if they are nonzero, have
opposite signs, retain every branch gate and exclude `s=0`.

Perform exactly 80 bisections.  At each step evaluate the midpoint, retain
the subinterval with opposite endpoint signs, and stop early only if `E` is
exactly zero at the working precision.  No Newton step, secant step, bracket
expansion, alternate seed or root library is permitted.  Use the final
midpoint as `s_root` and report the final interval width.

Require

```text
final width       < 2e-28,
abs(E(s_root))    < 1e-25,
abs(s_root)       > 1e-7,
s_root            < 0.
```

These gates establish one bracketed root but do not prove uniqueness inside
the bracket.

## 3. Branch and full-equation gates

At both endpoints, every midpoint and the final root require:

- all 2400 simplex representatives counted with one timelike Gram
  direction;
- positive leading-minor margin;
- minimum angle-argument modulus above `1e-6`;
- maximum action/gradient imaginary contamination below `1e-70`;
- positive diagonal, pole magnitude and upper boundary square.

At the final root require separately

```text
max abs(30 diagonal residuals) < 1e-60,
max abs(5 pole residuals)      < 1e-25,
max abs(all 35 residuals)      < 1e-25,
within-diagonal spread         < 1e-60,
within-pole spread             < 1e-60.
```

Failure of the diagonal gate is a local-variation refutation even if `E=0`.

## 4. Canonical-junction gate

At the same root compare all 30 newly evaluated pre-momenta with the stored
forward target.  Define

```text
r_p = p_pre(root)-p_target.
```

Require

```text
||r_p||_2 <= 10*cusp_uncertainty_norm,
max component spread of r_p <= 10*cusp_uncertainty_norm.
```

Also compare against the exact static formula
`+epsilon3*L0*tau/4` as a diagnostic, but use the independently stored target
for acceptance.  A uniform nonzero mismatch does not pass.

## 5. Parity and continuum diagnostics

Require the even and odd final root intervals to overlap after enlarging
each endpoint by `1e-28`; otherwise the outcome is schedule-dependent.
Also require

```text
abs(s_even-s_odd)                         < 1e-27,
abs(exp(s_even)-exp(s_odd))               < 1e-27,
max abs(p_pre_even-p_pre_odd)             < 1e-24,
max abs(p_post_even-p_post_odd)           < 1e-24.
```

Report, without using it as an acceptance target, the continuum control

```text
s_cont = -zeta^2*rho/(2 L0^2),
zeta=(pi^2*sqrt(2)/50)^(1/3),
```

and the ratio `s_root/s_cont`.  Agreement or disagreement with the continuum
does not rescue a failed Regge or canonical gate.

## 6. Mechanical outcome hierarchy

Assign exactly one outcome in this order:

1. `HOMOTHETIC_FORWARD_CONTROL_FAILED` if an upstream artifact, geometry or
   branch control fails;
2. `HOMOTHETIC_FORWARD_BRACKET_REFUTED` if either parity lacks the frozen
   sign change;
3. `HOMOTHETIC_FORWARD_ROOT_NUMERICALLY_OPEN` if bisection does not meet its
   width/residual gates;
4. `HOMOTHETIC_MINISUPERSPACE_ONLY` if `E=0` but any complete internal gate
   fails;
5. `HOMOTHETIC_STATIONARY_NOT_CANONICAL` if all 35 equations pass but the
   momentum junction fails;
6. `HOMOTHETIC_FORWARD_SCHEDULE_DEPENDENT` if both separate roots pass but
   their parity comparison fails;
7. `HOMOTHETIC_FORWARD_ROOT_ACCEPTED` only if every preceding gate passes.

## 7. Acceptance boundary and interpretation

Outcome 7 establishes a **DERIVED COMPUTATIONAL** non-static homogeneous
slab with:

- fixed conserved mass;
- all complete internal equations satisfied;
- contraction if `s_root<0`;
- correct canonical gluing to the published initial sandwich.

It is then legitimate to call this the first homogeneous forward tick of the
fixed carrier at the chosen lapse.  It is not legitimate to call the lapse,
time unit, branch, or full 65-variable evolution map unique.  Multi-tick
evolution, stability, continuum convergence, `c` and Planck scales remain
open.

Only the new targeted verifier will be run.  It must be registered exactly
once in `reproducible/run_all.py`; the full suite will not be run.
