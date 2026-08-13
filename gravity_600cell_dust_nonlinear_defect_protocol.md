# Preregistration: raw-residual nonlinear defect correction

Date: 2026-08-13

Prior-art gate: `81b1aa1`

First nonlinear result and post-result diagnostic: `b4a7828` plus the
committed amendment preceding this protocol

Status: **frozen before iterating a second defect-correction step or
validating any corrected candidate from this protocol**.

The preceding one-step reduction was inspected and is explicitly not blind
evidence.  This protocol tests whether the same fixed algorithm closes all
equations and reports all sixteen hit outcomes.

## 1. Fixed cases and unchanged scientific scope

Load without alteration the four boundary directions, two signs, two
schedules, amplitude `eta=1e-4`, final-boundary logarithmic vectors and lapse
grid

```text
t = 0,+0.05,+0.10,-0.05,-0.10
```

from protocol `80f8de7` and result `b4a7828`.  There are exactly sixteen
signed parity cases.  Do not recompute a favorable direction or reduce the
amplitude.

Use the same complete analytic logarithmic equation

```text
E=u*(partial S_total/partial u)/24,
F=Q^T E,
```

the same exact base complement `Q`, the same schedule-specific committed
quotient matrix `H_Q`, and the same parameterization `log u=log u_path(t)+Qz`.

## 2. Frozen defect iteration

For every grid point, start from the stored final `z` of the first nonlinear
run, not from a newly chosen initial vector.  At iteration `n`, compute

```text
p_n = H_Q^-1 F_n.
```

Try damping factors in the fixed order

```text
1,1/2,1/4,...,1/1024.
```

Accept the first trial `z-alpha*p_n` which is finite, remains on the certified
Lorentzian branch, and strictly reduces `norm(F)`.  The acceptance metric is
the physical transverse equation norm, not `norm(p_n)`; the post-result
diagnostic showed that the latter is unreliable under rotation of the four
soft directions.

Stop at the first of:

- `norm(F)<1e-12` and `norm(p_n)<1e-5`;
- no accepted damping;
- 50 accepted iterations.

Record every residual, correction norm, damping and branch margin.  A grid
point is `DEFECT_TRANSVERSE_SOLVED` only under the two numerical gates above.
No finite-difference Jacobian or other optimizer is allowed.

## 3. Scalar localization

Use the unchanged scalar

```text
g=w_path(t)^T E/norm(w_path(t))
```

and floor `1e-10`.  Candidate and bisection rules, interval width `1e-8`,
maximum 30 steps, arithmetic-mean midpoint initialization and `1e-6`
deduplication are identical to protocol `80f8de7`, except every transverse
solve uses the defect iteration above.

At a bisection midpoint initialize `z` by the mean of the two endpoint
solutions.  If either endpoint is unresolved, the interval is unresolved.

## 4. Independent validation

Use exactly the same complete-action 100-decimal derivative audit and steps

```text
2e-5,1e-5,5e-6
```

as protocol `80f8de7`.  Retain every gate unchanged:

- full equation, transverse correction and collective scalar are each within
  ten times their independently extrapolated errors;
- correction and correction error below `1e-5`;
- all 210 validation geometries pass branch gates;
- maximum imaginary contamination below `1e-80`;
- analytic/action equation agreement below `1e-8`.

No candidate passes merely because its raw binary64 transverse residual is
small.

## 5. Outcomes and evidence boundary

Use the same per-case outcomes and hit fractions over 16 signed cases, 8
direction-parity pairs and 4 contrasts.

Because a one-step improvement was inspected before preregistration, any
positive result has **PATTERN-INFORMED METHOD provenance**.  Robust evidence
requires a large hit fraction and action validation; one isolated hit remains
look-elsewhere limited.

If all sixteen scans resolve and none validates, the relative-phase sector is
killed at `eta=1e-4, |t|<=0.1` on this carrier.  If at least one scan remains
unresolved, the verdict is numerical OPEN.  No outcome establishes amplitude
scaling, a second slab, the full 840-edge carrier, a clock, speed limit or
Planck scale.
