# Control correction: adversarial refined boundary cotangent

Date: 2026-08-21

First-failure artifact and diagnosis are preserved before this correction.
Evaluate the kernel-control boolean inside `mp.workdps(100)`:

```text
abs(delta pullback)      < 1e-90,
abs(distance-1e-6)       < 1e-90.
```

Reuse that boolean outside the context.  Everything else remains frozen.
