# Preregistration: nonlinear continuation of the four relative-phase modes

Date: 2026-08-13

Prior-art gate: `81b1aa1`

Weak-scale path result: `ae902a8`

Precision linear response: `29a779f`

Status: **frozen before evaluating any deformed-boundary residual or running
any nonlinear solve**.

## 1. Complete question and fixed hypotheses

On the same order-24 invariant 600-cell carrier, fixed old regular boundary,
unrounded source `l0`, source dust mass and two derived even/odd schedules,
ask whether the four soft relative-phase directions of the linear response
integrate to complete nonlinear stationary solutions for a nonregular final
boundary.

The final-boundary perturbations preserve the geometric mean of the thirty
squared edge variables exactly in logarithmic coordinates.  The collective
lapse is not discarded as gauge after deformation: all 35 internal equations
must pass an independent action-only audit.

This is a local existence probe at one frozen predicted internal displacement
scale, not a proof for all 29 boundary directions or finite cosmological
evolution.

## 2. Prior-art boundary

Free lapse on a symmetric Regge background and lapse-dependent nonlinear
pseudo-constraints are **KNOWN STRUCTURE**.  Data free at linear order can be
fixed a posteriori by nonlinear consistency conditions.  The method below is
a finite-dimensional Lyapunov--Schmidt continuation, not new formalism.

What is **OPEN** is whether any or all of the four explicit relative-phase
responses on this dust-filled 600-cell carrier possess a nearby complete
nonlinear solution.

## 3. Frozen direction enumeration: N = 4

Do not select a largest singular vector.  The leading four response singular
values differ by only about `4e-6` relatively, so an individual SVD vector is
not a canonical direction.

Use the even-schedule precision result only to define the common physical
boundary directions tested in both schedules.  Let

```text
Y : 30 x 29 exact zero-sum boundary basis,
R : 35 x 29 internal linear response,
P5 = I5-ones(5)ones(5)^T/5,
L = P5 R_pole : R^29 -> ones(5)^perp.
```

The preregistered singular values of `L` are approximately

```text
4.742934e5, 4.742914e5, 4.742908e5, 4.742897e5, 6.98e-11,
```

and its absolute-`1e-8` rank is four.  Define the four ordered Helmert phase
contrasts

```text
c_k = (1,...,1,-k,0,...,0)/sqrt(k(k+1)),  k=1,2,3,4.
```

For each `k`, set

```text
a_k = pinv(L,rcond=1e-12)c_k,
y_k = Y a_k / norm(Y a_k).
```

These four directions are fixed by the derived five-stage order and the
minimum-norm inverse.  Test all four; no direction is dropped after seeing a
result.  Report the full vectors, their zero sums, the rank and all singular
values before any solve.

## 4. Frozen amplitudes and look-elsewhere count

For each direction let

```text
s_k = norm(R Y^T y_k).
```

The preregistered values lie between `4.7428998e5` and `4.7429293e5`.
Use exactly one predicted internal displacement scale

```text
eta = 1e-4,
boundary amplitude alpha_k = eta/s_k.
```

Thus `alpha_k` is about `2.1084e-10`.  Test both signs, all four directions
and both schedules:

```text
N_cases = 4 x 2 x 2 = 16.
```

The final squared boundary variables are

```text
v = l0^2 exp(sign*alpha_k*y_k).
```

No amplitude is changed after convergence behavior is known.  A later
scaling test at other `eta` values requires a new protocol.

## 5. Coordinates and transverse solve

At collective coordinate `t`, parameterize positive internal variables by

```text
rho(t)=tau^2 exp(t), q(t)=l0^2-rho(t),
log u(t,z)=log(q repeated 30,rho repeated 5)+Qz,
```

where `Q` is the committed exact base-tangent complement.  Let

```text
E=u*(partial S_total/partial u)/24,
F=Q^T E,
p=H_Q^-1 F.
```

Use SciPy `least_squares` on `p`, with the frozen options

```text
method='trf', jac='3-point', diff_step=1e-4,
xtol=ftol=gtol=1e-12, max_nfev=800, x_scale=1.
```

At `t=0`, initialize `z` with the committed linear coefficient response to
the signed boundary vector.  Scan in the order `0,+0.05,+0.10,-0.05,-0.10`,
using the nearest resolved solution as the next initial value on each side.

A transverse solve is resolved only if

```text
norm(F)<1e-9, norm(p)<1e-5,
```

all iterates evaluated by the residual remain finite, and the final geometry
is Lorentzian with Gram modulus above `1e-8` and angle-argument modulus above
`1e-6`.  Optimizer success without these physical residual gates does not
count.

## 6. Frozen scalar localization

At every resolved transverse point compute

```text
g=w_path(t)^T E/norm(w_path(t))
```

from the analytic complete gradient only for localization.  Use scalar floor
`1e-10` because this is not the scientific validation.

Candidates are every grid point with `abs(g)<=1e-10` and every adjacent
resolved pair with opposite signs and endpoint magnitudes above the floor.
Refine sign-changing intervals by bisection to width `1e-8` or
`abs(g)<=1e-10`, at most 30 steps.  Initialize each midpoint with the
arithmetic mean of endpoint `z` vectors and use the identical transverse
solver.  Deduplicate candidates separated by less than `1e-6` in `t`.

If a scalar is within the floor but the full action-only audit fails, it is a
false binary64 candidate, not a root.

## 7. Independent candidate validation at 100 decimals

For every localized candidate, reconstruct all 35 logarithmic derivatives
from the complete action at 100 decimal digits and steps

```text
h=2e-5,1e-5,5e-6
```

using the sixth-order construction in protocol `8380f0d`.  Let `E6` be the
row and `dE=E6-R23`.  Define

```text
delta=norm(H_Q^-1 Q^T E6),
epsilon_delta=norm(H_Q^-1 Q^T dE),
g6=w_path(t)^T E6/norm(w_path(t)),
epsilon_g=abs(w_path(t)^T dE/norm(w_path(t))).
```

A candidate is `NONLINEAR_STATIONARY_CANDIDATE` only if:

- `norm(E6)<=10 max(norm(dE),1e-30)`;
- `delta<=10 max(epsilon_delta,1e-30)`;
- `abs(g6)<=10 max(epsilon_g,1e-30)`;
- `epsilon_delta<1e-5` and `delta<1e-5`;
- all 210 validation geometries pass the Lorentzian branch gates;
- maximum imaginary contamination is below `1e-80`;
- binary64 and action-only equation rows agree in norm below `1e-8`.

This is a computational weak-scale certificate, not an interval proof of an
exact algebraic root.

## 8. Frozen outcomes and hit fractions

Per case assign:

- `NONLINEAR_CONTINUATION_FOUND_IN_SCAN` if at least one candidate passes;
- `NO_NONLINEAR_CONTINUATION_IN_FROZEN_SCAN` only if all five transverse grid
  solves and all sign-changing intervals resolve and no candidate passes;
- `NONLINEAR_CONTINUATION_NUMERICALLY_UNRESOLVED` otherwise.

Report all passing candidates; do not select only one lapse.  Report hit
fractions over:

```text
16 signed parity cases,
8 direction-parity pairs (both signs separately visible),
4 phase contrasts.
```

The route advances if at least one preregistered signed case has a validated
nonlinear continuation.  It establishes propagation only for that case.  A
failure in all sixteen resolved cases kills nonlinear continuation of the
complete relative-phase sector at `eta=1e-4` within `|t|<=0.1`, not outside
that scale or scan.

## 9. Claim boundary

No outcome proves continuum evolution, a graviton, multi-slab stability,
hyperbolicity, a speed limit, Planck units or coverage of the full 840-edge
carrier.  A positive result must next pass amplitude scaling and a second
slab.  A large hit fraction is evidence of a robust local solution map; one
isolated hit is **PATTERN / look-elsewhere limited**.
