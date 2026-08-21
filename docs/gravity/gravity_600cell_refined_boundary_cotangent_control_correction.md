# Control correction: refined boundary cotangent verifier

Date: 2026-08-21

First-failure commit: `5ed87f6`.

The first run's sole failure was a mixed-precision equality in the synthetic
kernel control.  Before rerunning, compute the control boolean inside a
100-decimal context and require

```text
abs(kernel_pullback-p_pre_fine) < 1e-90,
abs(kernel_distance-1e-6)       < 1e-90.
```

Both target numbers must be created inside that same context.  This changes
only how the already frozen negative control is represented.  The action,
all 24 schedules, `1e-60` scientific gates, finite-difference steps, raw and
normalized comparisons and outcome hierarchy remain unchanged.
