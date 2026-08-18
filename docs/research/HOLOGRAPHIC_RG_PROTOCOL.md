# Holographic RG Protocol

Objective: test whether the 600-cell boundary data support more than the
obvious `S^3` interpretation, namely a genuinely emergent `4D` bulk signal.

## Principle

Do not ask whether the model is "inspired by AdS/CFT".
Ask whether the spectral flow of the boundary data behaves as if a radial bulk
direction is present.

The first gate is intentionally severe:

- if the first gate fails, stop the gravity program in its holographic form
- if it passes, continue to universality and collective-mode tests

## Axiom 1 as a go/no-go test

Axiom 1:
the boundary spectral data must contain enough structured scale dependence to
support an emergent bulk interpretation.

Operational version:
the spectral flow must show a stable `4D-like` window, not only a `3D`
boundary window.

## Observable Set

### 1. Counting-function dimension

Use the positive spectrum of `D^2` and fit

`N(lambda) ~ lambda^(d/2)`

on rolling windows.

Interpretation:

- `d ~ 3` means the data behave like a boundary spectrum
- `d ~ 4` is the first sign that a bulk interpretation might be viable

### 2. Heat-flow effective dimension

Define

`K(t) = sum_{lambda>0} exp(-t lambda)`

and

`d_eff(t) = -2 d log K / d log t`

Search for stable windows where `d_eff` is approximately constant.

Interpretation:

- stable plateau near `3` -> boundary only
- stable plateau near `4` -> continue the holographic program

### 3. Coarse-graining universality

Only after step 1 is positive.

Apply 2-3 symmetry-respecting coarse-graining schemes and check whether the
dominant fitted coefficients stay approximately invariant.

If they drift strongly, the bulk interpretation is scheme-dependent and should
be rejected.

### 4. Collective-mode separation

Only after steps 1 and 2 are positive.

Perturb the boundary operator and compute the Hessian of the spectral action.
Look for a small set of low modes separated from the microscopic bulk of the
spectrum.

Without such a gap, there is no credible route to emergent geometry.

## Stop / Continue Criteria

`CONTINUE` only if both of the following hold:

- counting-function windows support `d = 4` within tolerance
- heat-flow windows support `d = 4` within tolerance

`STOP` if both of the following hold:

- counting-function windows support `d = 3`
- heat-flow windows support `d = 3`

`INCONCLUSIVE` for any mixed result.

## Current Command

Run:

```powershell
python reproducible\test_holographic_rg.py
```

This script is exploratory by design. It is allowed to say `STOP`.

That is a useful result, not a failure of the code.
