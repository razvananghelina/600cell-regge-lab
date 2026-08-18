# Nonlinear relative-phase continuation: first frozen result

Date: 2026-08-13

Prior-art gate: `81b1aa1`

Frozen protocol: `80f8de7`

Implementation: `12da490`

Registered verifier:
`reproducible/verify_gravity_600cell_dust_nonlinear_relative.py`

Machine-readable result:
`reproducible/gravity_600cell_dust_nonlinear_relative.json`

Targeted run: **10/10 implementation checks passed**.  The full suite was not
run.

## 1. Frozen verdict

```text
OPEN NUMERICALLY.
```

The hit fractions are:

```text
validated signed cases       0/16,
direction-parity pairs       0/8,
phase contrasts              0/4.
```

This is not a kill.  All sixteen cases received

```text
NONLINEAR_CONTINUATION_NUMERICALLY_UNRESOLVED
```

because none completed the entire frozen transverse scan.  Therefore the
protocol did not authorize `NO_NONLINEAR_CONTINUATION_IN_FROZEN_SCAN`.

## 2. Direction enumeration was nonarbitrary

The map from the 29 zero-sum boundary directions to the relative five-pole
sector has singular values

```text
474293.4383, 474291.4080, 474290.7555, 474289.6878, 6.98e-11
```

and absolute-`1e-8` rank four.  The first four values are nearly degenerate,
confirming that choosing one largest-SVD vector would have been arbitrary.
All four ordered Helmert phase contrasts and both signs were tested in both
schedules.

The boundary logarithmic amplitudes were approximately `2.1084e-10`, fixed so
that the predicted internal displacement norm was `eta=1e-4` in every case.

## 3. What happened in the even schedule

At `t=0`, all eight even cases passed the binary64 transverse localization
gates, producing one candidate each.  The other four lapse-grid points per
case did not pass the absolute transverse `1e-9` gate, although every final
geometry remained Lorentzian.

All eight `t=0` candidates were then rejected by the independent
100-decimal action audit.  Representative ranges were:

```text
action equation norm         9.69e-11 ... 4.36e-10
action empirical error       about 2.64e-14
quotient correction proxy    4.45e-8  ... 9.71e-8
correction-proxy error       about 1.71e-12
collective scalar            about 9.8e-17
collective scalar error      about 8.5e-19
binary/action row difference 3.8e-13 ... 5.6e-13
```

Thus the rejection is not caused by disagreement between the analytic
gradient and action-only implementation.  The residuals are resolved nonzero
under the frozen action estimator.  Yet the remaining predicted correction
is below `1e-7`, only about `0.1%` of the imposed internal displacement.  This
is evidence that the frozen finite-difference least-squares solver stopped
prematurely near a possible solution, not evidence that no solution exists.

## 4. What happened in the odd schedule

All forty odd transverse grid solves remained above the frozen absolute
`1e-9` residual gate, despite optimizer success flags and valid Lorentzian
geometries.  Therefore no scalar candidate was localized and no odd
action-only validation was triggered.

Typical odd transverse residual norms were `1e-8 ... 2e-7`, with
preconditioned correction norms mostly `1e-6 ... 2e-5`.  This is the same
binary64 weak-mode sensitivity seen in earlier Hessian work, now in a
nonlinear solve.

## 5. Interpretation

- **NEGATIVE under frozen solver:** no candidate passed the complete action
  gates.
- **OPEN physically:** existence or nonexistence of a nonlinear continuation.
- **NOT A KILL:** every signed case had an incomplete frozen scan.
- **DERIVED:** the simple `3-point` SciPy solve is not accurate enough for
  these `4.6e-8` soft modes at `eta=1e-4`.
- **NOT SUPPORTED:** any claim that a shape perturbation already constitutes
  a physical tick or evolution.

## 6. Correct next correction

Do not loosen the action gates and do not reduce `eta` after seeing the
failure.  Replace only the localization algorithm:

1. use the complete analytic equation row, whose agreement with the
   action-only row is now measured below `6e-13`;
2. take explicit Newton/defect-correction steps with the committed quotient
   Hessian rather than allowing a numerical `3-point` Jacobian to declare
   optimizer convergence;
3. accept steps only when the preconditioned residual decreases;
4. retain the same 16 cases, `eta=1e-4`, lapse grid, scalar rule and
   100-decimal candidate validation.

A post-result one-step diagnostic may decide whether this correction is
promising, but it cannot count as evidence until separately preregistered.

## 7. Post-result one-step diagnostic

The non-evidential diagnostic was run on all sixteen `t=0` grid states.  One
explicit defect step

```text
z_new = z-H_Q^-1 Q^T E
```

reduced the physical transverse equation norm as follows:

```text
even: 9.72e-11 ... 4.36e-10  -> 6.48e-13 ... 9.27e-13,
odd : 1.21e-8  ... 1.39e-7   -> 3.03e-11 ... 8.00e-11.
```

However, the norm of the next fixed-Hessian correction did not decrease
monotonically in every case.  This distinguishes two metrics:

- `norm(Q^T E)` is the physical equation residual and fell by orders of
  magnitude;
- `norm(H_Q^-1 Q^T E)` is highly sensitive to small rotation/error in the
  four soft eigenvectors and can grow even when the equations improve.

This diagnostic supports a separately preregistered raw-residual line-search
correction.  It is **PATTERN / method diagnostic**, not a nonlinear hit.
