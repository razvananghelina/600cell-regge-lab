# Preregistration: stationary-base Lyapunov--Schmidt audit

Date: 2026-08-13

Prior-art and framing correction: `86f6ce7`

Precision quotient result: `29a779f`

Published dust control: `66a6465`

Status: **frozen before solving any transverse equation away from the
printed base or inspecting the scalar reduced equation**.

## 1. Question and complete hypotheses

For each of the two already derived even/odd five-stage schedules, keep fixed:

- the old and final regular 600-cell boundaries;
- the published unrounded `l0` reconstructed from its displayed formula;
- the published dust mass reconstructed from its displayed formula;
- the printed reference `tau=0.0102` used only to parameterize the collective
  path;
- the same 35 orbit-reduced internal carrier and real Lorentzian branches.

Does the complete 35-equation dust Regge system possess a genuinely
stationary point near the action-flat collective path, once the four soft
relative-pole directions are resolved on their own `4.605e-8` scale?

This protocol does not perturb a boundary, compare to a physical target,
fit the dust mass, or assume that collective lapse remains gauge.

## 2. Known mechanism and open object

Higher-order consistency conditions that fix background lapse parameters and
turn linear constraints into pseudo-constraints are **KNOWN** in Regge
calculus.  The method below is a finite-dimensional Lyapunov--Schmidt
reduction, not a new formalism.

The existence, location and multiplicity of stationary points for this
specific dust-filled, five-stage 600-cell carrier are **OPEN**.  In
particular, the earlier `1e-7` absolute stationarity gate was not conditioned
on the four much smaller transverse curvatures and cannot answer this
question.

## 3. Frozen coordinates

For collective coordinate `t`, use the action-flat path

```text
rho(t) = tau^2 exp(t),
q(t)   = l0^2-rho(t).
```

Let `u_path(t)` be the 35-vector of logarithms of thirty copies of `q(t)` and
five copies of `rho(t)`.  Load the exact normalized base tangent and its
deterministic Householder complement `Q` from protocol `da34272`.  Parameterize
all nearby positive internal variables by

```text
log u(t,z) = u_path(t) + Q z,   z in R^34.
```

The path tangent at fixed `z` is

```text
w_path(t) = (-rho(t)/q(t) repeated 30, 1 repeated 5).
```

For the complete analytic action gradient define the per-edge logarithmic
equations

```text
E(u) = u * partial S_total/partial u / 24.
```

The reduced system is

```text
F(t,z) = Q^T E(u(t,z)) = 0          (34 equations),
g(t)   = w_path(t)^T E(u(t,z(t)))
         / norm(w_path(t)) = 0      (one equation).
```

Because `[Q,w_path(t)]` is invertible throughout the frozen interval, both
conditions together are equivalent to all 35 local equations.  The scalar
equation is not discarded as gauge.

## 4. Frozen transverse solver and scan

Use the precision-corrected base quotient matrix `H_Q` from commit `29a779f`
only as a fixed preconditioner.  At fixed `t`, iterate

```text
delta z = -H_Q^-1 F(t,z).
```

Accept the first damping factor in

```text
1, 1/2, 1/4, ..., 1/1024
```

which strictly lowers `norm(H_Q^-1 F)`.  Stop after at most 100 accepted
iterations.  A transverse solve is resolved only if

```text
norm(F) < 1e-12
and norm(H_Q^-1 F) < 1e-5.
```

No alternative optimizer or retuned damping sequence is allowed after seeing
the result.  Failure is `TRANSVERSE_SOLVE_UNRESOLVED`.

Evaluate the fixed grid

```text
t = -0.10, -0.075, -0.05, -0.03, -0.02, -0.01,
     0,
     0.01, 0.02, 0.03, 0.05, 0.075, 0.10.
```

Start at `t=0,z=0`; continue separately toward positive and negative `t`,
always using the nearest resolved `z` as the next initial point.  Report every
`z`, residual, damping history, branch margin and scalar `g`; no grid point is
discarded.

## 5. Frozen scalar-root rule

The scalar resolution floor for the binary64 localization stage is `1e-12`.
Candidate roots are:

1. every resolved grid point with `abs(g)<=1e-12`;
2. every adjacent resolved pair whose scalar values have opposite signs and
   both exceed that floor in magnitude.

Refine every sign-changing interval by bisection, not interpolation or a
derivative fit.  At each midpoint solve the transverse equations by the same
frozen iteration, initialized at the arithmetic mean of the two endpoint
`z` vectors.  Replace the endpoint having the midpoint's sign and continue.
Stop at interval width `1e-10` or `abs(g)<=1e-12`, with a maximum of 40
bisections.  Deduplicate candidates separated by less than `1e-8`.  Do not
choose only the root nearest the printed `tau`; multiplicity is a scientific
result.

All trial geometries must remain Lorentzian with minimum Gram modulus above
`1e-8` and minimum angle-argument modulus above `1e-6`.  A branch failure is
reported and makes the affected interval unresolved.

## 6. Independent 100-decimal action-only audit

For every localized candidate, reconstruct the 35 logarithmic derivatives
directly from the complete action at 100 decimal digits.  At each of the
frozen steps

```text
h1=2e-4, h2=1e-4, h3=5e-5
```

and each internal coordinate `i`, evaluate

```text
D_h[i] = [S(log u+h e_i)-S(log u-h e_i)]/(48 h).
```

Form the same fourth- and sixth-order extrapolations

```text
R12=(4 D_h2-D_h1)/3,
R23=(4 D_h3-D_h2)/3,
E6 =(16 R23-R12)/15,
epsilon_action=norm(E6-R23).
```

Report the complete row, all imaginary parts and branch margins.  Define

```text
action_floor=max(epsilon_action,1e-30),
delta_transverse=norm(H_Q^-1 Q^T E6).
```

A candidate is `WEAK_SCALE_STATIONARY` only if all hold:

- `norm(E6) < 1e-10`;
- `abs(w_path^T E6/norm(w_path)) < 1e-11`;
- `delta_transverse < 1e-5`;
- the sixth-order imaginary contamination is below `1e-80`;
- all 210 action perturbations remain on the certified branches;
- the binary64 candidate and action-only `E6` agree to norm below `1e-9`.

These gates certify a base point to a scale adequate for the next proposed
`1e-4` internal perturbations.  They are not an interval proof of an exact
root.  Also report `norm(E6)/action_floor`; do not require it to resemble one.

## 7. Frozen outcomes

Assign one outcome per schedule:

- `UNIQUE_WEAK_SCALE_STATIONARY_BASE_IN_SCAN`;
- `MULTIPLE_WEAK_SCALE_STATIONARY_BASES_IN_SCAN`;
- `NO_STATIONARY_BASE_IN_FROZEN_SCAN` if all grid/interval solves are
  resolved and no candidate passes;
- `STATIONARY_BASE_NUMERICALLY_UNRESOLVED` otherwise.

The two parities need not agree.  Root count outside `[-0.1,0.1]` is not
claimed.

## 8. Claim boundary and decision rule

A resolved stationary base would show whether the printed `tau` is merely a
nearby coordinate choice or is shifted/fixed by a pseudo-constraint.  It
would not yet validate the previously computed 29-direction response: the
Hessian, exact-null status and boundary mixed row must be recomputed at the
new base before any nonlinear boundary continuation.

If no stationary base exists in the frozen scan, the published sandwich does
not provide a local dynamical base on this restricted carrier.  If the audit
is numerically unresolved, no physical conclusion is drawn and no threshold
is loosened.  No outcome selects a clock, Planck scale, graviton or continuum
limit.
