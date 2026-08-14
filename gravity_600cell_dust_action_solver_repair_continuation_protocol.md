# Preregistration addendum: calibrated complete-action solver continuation

Date: 2026-08-14

Prior-art gate: `6b7f9e4`

Base solver protocol: `4b6b10c`

Failed calibration result and post-result search: `87c289f`

Status: **frozen before evaluating any deformed target row**.  The failed run
evaluated exactly zero target rows.

## 1. Exact amendment

Protocol `4b6b10c` is retained verbatim except for section 4, calibration gate
3.  Replace

```text
maximum absolute agreement with the old 60-decimal control < 1e-10
```

by

```text
maximum absolute agreement with the old 60-decimal control < 5e-8.
```

The replacement value is not the observed discrepancy `3.544337829e-9`.
It is the threshold already registered by the upstream control for its own
60-decimal `3e-6` action-only difference versus the analytic gradient.  The
upstream measured errors were `1.106e-8` and `1.262e-8`, so asking another
calculation to agree with that reference at `1e-10` was not a valid control.

The old failed result remains immutable and reproducible.  This continuation
writes a distinct checkpoint, JSON artifact and verifier name.

## 2. Everything that remains frozen

No other mathematical or numerical choice changes:

- the same complete Regge-plus-dust action and `35=30+5` carrier;
- the same 80 `(parity,direction,sign,t,boundary)` states;
- the same `final_z` starts from `64a13f6`;
- 100-decimal variables, actions, differences and projections;
- operational steps `(1e-20,1e-15)`;
- disjoint validation steps `(3e-20,3e-15)`;
- the measured shadow-difference self-agreement gate;
- all branch and `1e-70` imaginary gates;
- the binary Jacobian only as proposal/preconditioner;
- its steps, singular-value, step-change and model-error gates;
- the active natural/raw merit switch;
- damping `1,...,1/1024` and twelve accepted iterations;
- every transverse-zero, validation and scalar threshold;
- the five-point grid, forced brackets and `16/8/4` hit accounting;
- every acceptance, kill and claim boundary in `4b6b10c`.

In particular, neither the positive scalar pattern nor a target residual was
used to choose this amendment.

## 3. Calibration and outcome rule

The symmetric control is again evaluated before loading or evaluating a target
state.  If either parity fails any gate, including the inherited internal
self-agreement gate or the amended upstream-accuracy gate, stop with
`DERIVATIVE_CALIBRATION_FAILED` and zero target rows.

If calibration passes, execute the target solver exactly as already frozen.
The scientific outcome hierarchy is unchanged.

## 4. Evidence boundary

Passing calibration is a **DERIVED CONTROL**, not evidence for a physical
tick.  A validated target root and reduced scalar retain the labels of
`4b6b10c`.  Any unresolved state keeps the route **OPEN NUMERICALLY**.

The amendment repairs a demonstrably invalid comparison tolerance.  It does
not establish physical time, dynamics beyond one restricted slab, inertia,
mass, a causal speed limit, Planck units or particle masses.

