# Second-run correction: missing sector-helper constant

Date: 2026-08-19

After the sparse-import repair, the second adversarial execution passed
provenance, the primary-result control, the `43/43` geometry import and the
independent global ranks `470/354/4`.  It then stopped before sector or target
construction because the AST-extracted `high_precision_sector_bases` helper
refers to the source module's global imaginary unit `I`:

```text
NameError: name 'I' is not defined
```

No adversarial residual or JSON artifact was produced.  The frozen repair is
only

```text
I = mp.mpc(0, 1).
```

No numerical method, carrier, matrix, threshold or outcome rule changes.
