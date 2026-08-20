# Preregistered correction: Lorentzian branch and reality gate

Date: 2026-08-20

Second failure commit: `72c1749`.

The second execution proved that the original generic imaginary-curvature
gate was not a legitimate condition in the complex Lorentzian Regge
convention. At the corrected anchors the complete action and all ten equations
are real to `1e-77`, the angle identities and logarithm arguments are healthy,
but individual hinge curvatures have expected imaginary components.

## Frozen correction

Define

```text
maximum_physical_imaginary = max(
    |Im complete action|,
    |Im gravitational action|,
    max_i |Im G_i|
).
```

Use this value, not `maximum_imaginary_curvature`, in every trial-point,
anchor, refinement and final-root reality gate. Continue to report
`maximum_imaginary_curvature` separately as a branch diagnostic.

At 50-digit solver precision require

```text
minimum angle-argument modulus > 1e-10,
maximum angle-identity residual < 1e-20,
maximum_physical_imaginary < 1e-18.
```

At the 140-digit final validation replace the protocol line

```text
max imaginary action/gradient/curvature < 1e-50
```

by

```text
max imaginary action/gradient < 1e-50.
```

This correction restores the already certified Lorentzian convention; it does
not relax reality of the action or equations. The corrected anchors from
commit `69ca0f3` remain unchanged. All 120 attempts, bounds, seeds, solver
options, root residual gates and outcomes remain unchanged.

If the action or equations acquire a resolved imaginary part, or an angle
argument approaches zero, the point still fails the branch gate.
