# Pre-output precision correction for nonlinear covariance seeds

Date: 2026-08-17

## What happened

The Stage-A implementation was committed as `fc45d7e` before its first run.
During that run, after the even Hessian but before the odd Hessian completed,
the console reported

```text
s_min(J)          = 4.24456181e-9
epsilon_J         = 1.08689e-16
response_error    = 4.15037e-16.
```

The rank remained decisively resolved, but the calibration was seven orders
weaker than the independently certified dynamic tangent calculation, whose
`epsilon_J` is about `6.8e-23`.

Inspection found that `DERIVATIVE_STEPS`, `ETA` and the arithmetic floor were
constructed as `mpmath.mpf` values before setting `mp.dps=100`.  They therefore
carried the default-context initialization into the high-precision run.  The
audited tangent verifier sets the precision first.

The process was interrupted during the odd derivative batch.  It exited 130,
wrote no seed artifact, derived no amplitude or case list, evaluated no
nonlinear perturbed case and compared no nonlinear output.

## Frozen correction

Make exactly one semantic change before rerunning Stage A:

```text
set arb.mp.dps = DPS immediately after defining DPS,
then construct every arb.mpf protocol constant.
```

Do not alter input hashes, derivative values, gates, directions, displacement
scale, case count, amplitudes, classifiers or output fields.  Commit the
corrected implementation before rerunning.

## Status

- **DERIVED IMPLEMENTATION FAILURE:** the interrupted run did not reproduce the
  intended arbitrary-precision initialization order.
- **CONTROL:** no evidential nonlinear result or candidate artifact exists from
  that run.
- **OPEN:** the corrected Stage-A seed enumeration.
