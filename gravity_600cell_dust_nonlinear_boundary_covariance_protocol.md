# Preregistration: nonlinear anisotropic boundary covariance

Date: 2026-08-17

Prior-art gate: `526a202`.

Status: frozen before recomputing the dynamic response matrix, deriving a new
amplitude, or solving a perturbed nonlinear canonical equation.

## 1. Frozen inputs

Require the following SHA-256 values:

```text
accepted dynamic tick
4b1c59c0518eec11b88b140cdecdf558d762c0d70b4826a758f67544e14ac5b9

blind dynamic tangent artifact
1ed8d63b4c8a6a4530570a2894820962c7c3c7852747a1112cdf1b242253dbb5

target-blind boundary candidate enumeration
51b52457eba84ca1e41926b6e4fb1c51032f788b70bde916a3fb755d0323cb3e

older nonlinear relative-direction artifact
6e7d108ec7b1a2c80b412134a301084aea14f9457fedba1fd840820ad6f558dd

audited canonical-Hessian source
396c491fe51a9f5e04fa8402e2e5b16884fe23fc5057d8ded325e6064fbd3b9e

one-slab carrier/action source
ad93cdd08fabeeee56b009f23936696837c4362f88ae23f92a36d0395e61ffaf.
```

Require the four older directions, their rank-four certificate and the
previously frozen displacement scale

```text
ETA = 1e-4.
```

Use the unique boundary candidate tagged
`IDENTICAL_PHYSICAL_EDGE_SETS`; no other one of the 60 passing permutations may
be selected.

No continuum spectrum, desired nonlinear answer, speed, experimental number or
full-carrier result may be loaded.

## 2. Stage A: calibrated response and frozen case list

At `DPS=100`, reconstruct each accepted dynamic background and the complete 95
by 95 logarithmic Hessian using exactly the already certified derivative pairs:

```text
operational primary = 1e-20,
operational shadow  = 1e-15,
validation primary  = 3e-20,
validation shadow   = 3e-15.
```

Retain all old controls: every displaced action evaluation must stay on the
Lorentzian branch; every operational-primary versus validation-primary Hessian
entry must pass the calibrated gate; Hessian reciprocity must pass; and all 65
singular values of

```text
J = [[ K_XX,  K_XN],
     [-K_OX, -K_ON]]
```

must be resolved nonzero.

For canonical input `a=(delta o,delta p_pre)`, use

```text
R = [[-K_XO, 0],
     [ K_OO, I]],

Y = J^-1 R,
```

so `Y a=(delta x,delta n)` is the linear unknown response.  Derive both
operational and validation response matrices independently.

Read each stored 30-vector direction, subtract its mean and renormalize it at
100 decimals.  Let

```text
p_star = abs(mean(p_pre_even_base)).
```

For direction `d`, define two dimensionless unit input rays:

```text
POSITION: (d,0),
MOMENTUM: (0,p_star*d).
```

Map the even ray to odd coordinates only with the preregistered physical-edge
permutation.  For each direction and sector, compute the Euclidean norm of the
65-component linear unknown response for both parities and both operational and
validation matrices.  Set one common amplitude

```text
amplitude = ETA / maximum(response norms).
```

Freeze exactly

```text
4 directions x 2 sectors x 2 signs x 2 levels = 32 paired cases,
levels = {1/2,1}.
```

Stage A must write every input, amplitude, response seed and control to a JSON
artifact with

```text
nonlinear_perturbed_action_evaluations = 0,
nonlinear_outputs_compared = false.
```

Commit that artifact before implementing or running Stage B.

## 3. Stage B: complete nonlinear canonical solves

For each frozen paired case and parity, set

```text
log q_old = base + signed_level*amplitude*delta_o,
p_target  = base + signed_level*amplitude*delta_p,
```

and seed the 65 logarithmic unknowns `(log x,log q_new)` with the corresponding
linear response.  Solve the complete real parts of

```text
gradient_internal[35] = 0,
-gradient_old[30] - p_target[30] = 0
```

using fixed-J Newton with the separately stored operational and validation
base matrices.  Use at most 20 iterations and backtracking factors
`1,1/2,...,2^-12`; accept a trial only if it stays on the certified Lorentzian
branch and satisfies the Armijo decrease

```text
||r_trial||_infinity <= (1-alpha/4)||r||_infinity.
```

Require for both independent solves:

```text
residual infinity norm < 1e-55,
maximum imaginary contamination < 1e-70,
one negative Gram direction in every simplex,
minimum leading minor > 0,
minimum angle-argument modulus > 1e-6,
final fixed-J correction infinity norm < 1e-45.
```

No failed case may be silently retried with another algorithm, amplitude or
direction.

## 4. Dimensionless outputs and calibrated comparison

For every converged solve define

```text
F = (log q_new[30], p_post[30]/p_star).
```

Map the even output to odd coordinates with
`C(Q)=diag(Q,Q)` from the unique physical-edge identification.  For each case
use the operational covariance defect

```text
d = ||F_odd - C(Q)F_even||_2.
```

Define its empirical uncertainty before classification as

```text
u = ||F_even_operational-F_even_validation||_2
  + ||F_odd_operational-F_odd_validation||_2
  + correction_output_even
  + correction_output_odd
  + 1e-70,
```

where `correction_output` is the norm of the output change after one final
fixed-J correction and re-evaluation.

Classify each of all 32 paired cases:

- `COVARIANT` if `d <= 10*u`;
- `BROKEN` if `d > 100*u`;
- `OPEN` otherwise.

Report every defect, uncertainty and the exact hit fractions.  If both levels
of a direction/sector/sign are `BROKEN`, report only as a diagnostic

```text
observed_order = log2(d_full/d_half).
```

Label `QUADRATIC_COMPATIBLE` only when the propagated defect intervals exclude
zero and the order interval lies inside `[1.5,2.5]`.  This diagnoses the first
breaking order; it is not a physical scaling law.

## 5. Mechanical outcome

- any `BROKEN`:
  `NONLINEAR_BOUNDARY_COVARIANCE_BROKEN_ON_FROZEN_CASES`;
- none broken, at least one `OPEN` or failed solve:
  `NONLINEAR_BOUNDARY_COVARIANCE_OPEN`;
- all 32 `COVARIANT` with every control passing:
  `NONLINEAR_BOUNDARY_COVARIANCE_CONSISTENT_ON_FROZEN_CASES`;
- any provenance, Hessian, branch or implementation control failure:
  `NONLINEAR_BOUNDARY_COVARIANCE_CONTROL_FAILED`.

Passing all 32 is **PATTERN / finite computational evidence**, not a theorem on
an open neighbourhood.  One resolved broken case falsifies nonlinear
covariance for the present quotient.

## 6. Scope

This protocol covers four precommitted shape directions, pure position and pure
momentum rays, their signs and two finite amplitudes.  It does not cover the
other 25 shape directions, mixed rays or cross-quadratic polarizations.  It
also does not cover the full 720-edge carrier, matter perturbations, refinement,
gravitons, a causal cone or a physical tick duration.

Register each stage before its first execution.  Run only the two targeted
verifiers, never the full suite.
