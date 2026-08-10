# Refined Whitney/Hopf verdict

Date: 2026-08-10

## Provenance

The refinement rule and the numerical result both preceded this comparison:

- tensor-definition commit: `be5914a`;
- blind-result commit: `973b512`;
- blind JSON SHA-256:
  `18bb2956b1a9d9acbe46fde8355afbba2b38404ac7273cb142b5ddeaa1b6c6d2`.

The rule used one preregistered fiber edge per parent tetrahedron to define

`Q_f=(w/V)e_f tensor e_f`, `Q_c=P_tangent-Q_f`,

and integrated those same tensors over all 24 barycentric children.  No fine
edge labels or level coefficients were fitted.

## Result

```text
                         coarse              first refinement
fiber kernel                12                       12
cross kernel                 1                        1
fiber gap             0.5647100088             0.5221352015
cross gap             2.8235500443             2.7875963884
cross/fiber ratio      5.0000000005             5.3388401713
```

The ratio drifts upward by `6.7768%`.  The fiber gap falls by approximately
`7.54%`, whereas the cross gap falls by only approximately `1.27%`.  Both
first positive clusters retain multiplicity four, and the kernel counts stay
fixed.

The failure is not numerical leakage from the inductive construction:

- mass compression residual: `1e-15`;
- fiber-form compression residual: `1e-15`;
- cross-form compression residual: `1e-15`;
- `Q_f` has exact observed nonzero eigenvalue `1/2`;
- `Q_c` is positive on the tangent space with minimum eigenvalue `1/2` to
  numerical precision.

## Verdict

**DERIVED NEGATIVE:** the equality of the base-level first-gap ratio with
`a_1=5` is not stable when the uniquely selected local fiber direction is
extended to the new Whitney refinement modes.

The base equality remains an exact theorem, but its evidential status is now
clearer: it is a property of the coarse graph spectrum preserved by the
coarse Whitney mass.  It is not a refinement-invariant propagation constant.

Therefore, within this preregistered local-tensor extension, the claim

`c^2=a_1=5`

is **CLOSED/KILLED when `c^2` means the first cross/fiber kinetic-gap ratio of
this inductive Whitney dynamics**.  This is not a uniqueness theorem for all
possible continuum extensions of a discrete Hopf edge labelling.

This does not refute the arithmetic bootstrap `a_1=5`, the 600-cell carrier,
or the existence of the Whitney dynamics.  It refutes one proposed physical
interpretation linking them.

## Framing attack

Could another spectral ratio remain five?  Almost certainly many ratios or
Rayleigh quotients can be formed after seeing the refined spectrum.  Searching
them now would recreate the look-elsewhere problem.  A replacement observable
is admissible only if independently motivated and preregistered before its
comparison with the seed.

Could the ratio return to five at deeper levels?  One level cannot exclude a
later limit numerically.  But exact equality and level-independence -- the
stated mechanism for a fixed fundamental `c` -- have already failed.  A new
claim about an asymptotic limit would require a new convergence hypothesis and
at least a second level; it cannot rescue the original claim retroactively.

## Consequence for the theory

- **DERIVED:** `a_1=5` still selects the Fibonacci bootstrap output.
- **DERIVED WITH STATED AUXILIARY SELECTORS:** it still leads to the 600-cell
  realization.
- **DERIVED:** the Whitney Galerkin-inductive dynamics still exists.
- **DERIVED NEGATIVE:** the new dynamics does not carry a fixed first-gap
  speed ratio equal to five.
- **OPEN:** physical `c`, Lorentzian time and the temporal/spatial
  normalization must come from a different dynamical mechanism.
