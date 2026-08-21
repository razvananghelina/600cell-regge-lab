# Frozen correction: primary serialization comparison

Date: 2026-08-21.

Adversarial protocol commit: `e8628bc`.

Registered adversarial verifier: `34fc405`.

Preserved first adversarial artifact: `c1371a6`.

First adversarial artifact SHA-256:
`edf372d841c6f1d69dfc521244a6aa37a2471f4adb6243c410c94ed249e9167e`.

Status: frozen after the first `7/9` run and before modifying or rerunning the
adversarial verifier.

## Failure

All direct full-action, precision, seed, reverse and hostile checks passed.
Only the final primary comparison failed.  The verifier required absolute
differences below `1e-70`, but the primary artifact serialized each root with
`mp.nstr(value,60)`.  Its information ceiling is therefore approximately
`1e-59` for order-one values.  The observed differences were

```text
branch A: 4.09e-60, 1.53e-62, 1.96e-60,
branch B: 2.06e-62, 1.67e-59, 3.54e-61.
```

This comparison could not pass even for the same underlying high-precision
root.

## Allowed correction

Keep every action, state, seed, direct solve, residual threshold, precision
run, physical inequality and outcome unchanged.  For the final comparison
only, require

```text
mp.nstr(independent_value,60)==primary_serialized_string
```

for `(h2,q2,scale_ratio)` on both branches.  This uses every digit retained by
the frozen primary artifact.  Also continue to print the absolute difference;
it may not exceed the primary serialization scale.

Before freezing this correction, all six 60-digit strings were checked and
matched character for character.  This is a disclosed representation repair,
not new physical evidence.  Any failure after the repair leaves the result
`FINITE_HEIGHT_TWO_SLAB_OPEN`.
