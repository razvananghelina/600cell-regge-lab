# Constrained H4 auxiliary failures are resolved

Date: 2026-08-21

Status: clean diagnostic result; a corrected adversarial gate is licensed but
has not yet been executed.

Protocol commits: `af8aa91`, `9e44dd7`.  Verifier registration: `c991173`.

The verifier passed `13/13` twice and produced the byte-identical artifact

```text
reproducible/gravity_600cell_refined_h4_constrained_response_auxiliary_diagnostic.json
SHA-256 f66177326afc3b3457a60b544745b739cbaa6b6d6e7f367b57d60f31eeeddeb7
```

No full suite, root search or deferred nonlinear census was run.

## Stationarity failure

For all ten internal directions on all 24 schedules, at both 180 and 220
decimal digits, the direct scalar-action ladder changed the maximum residual
from

```text
legacy fourth-order value  5.1661660e-54
tenth-order final value     5.9097062e-77.
```

Every one of the 480 precision-specific zero gates passes.  The maximum used
zero-envelope fraction is `5.9097e-17`; the maximum 180/220 precision fraction
is `7.9027e-100`.  A known even action is accepted, while adding a resolved
linear derivative `1e-20` is recovered and rejected by the zero gate.

**DERIVED COMPUTATIONAL:** the first audit's `1e-60` stationarity failure was
the truncation error of its under-resolved fourth-order scalar derivative.  It
is not a resolved nonzero internal gradient.

## Individual-curvature reality failure

On all 24 schedules and every one of the 210 direct-audit directions, the
global individual-curvature maxima at successively halved steps are

```text
2.0750292e-8, 1.0375146e-8, 5.1875730e-9,
ratios 2.0000000, 2.0000000.
```

The background maximum is `3.4262e-178`.  All `5040/5040` finest-step
plus/minus pairs have resolved smooth odd-leading behaviour; the maximum
even/odd ratio is `1.6772e-10`.  Simultaneously,

```text
relative imaginary complete action  9.1262e-177
angle-identity residual              8.6673e-178
minimum angle argument               0.91222028.
```

The smooth complex control passes with halving ratios `(2,2)`, while a
discontinuous branch control gives `(1,1)` and is rejected.

**DERIVED COMPUTATIONAL / STRUCTURAL:** individual Lorentzian curvature terms
become complex smoothly under arbitrary off-shell real log-edge
perturbations, vanish back to the real static background, and combine into a
real complete action on a safe analytic branch.  Requiring each individual
off-shell term to remain real was an overstrong category error.  The relevant
branch gates are continuity/safe analytic arguments and reality of the
complete action, not termwise curvature reality away from shell.

## Verdict and limits

The frozen diagnostic outcome is

```text
REFINED_H4_CONSTRAINED_RESPONSE_AUXILIARY_FAILURES_RESOLVED.
```

This result does not retroactively turn the first adversarial execution into
a pass.  It licenses a separately preregistered corrected gate that preserves
the already frozen direct matrices and replaces only the two falsified
auxiliary controls with the diagnostic evidence above.

The direct single-class result remains **OPEN** until that corrected gate
passes.  Even then it concerns only the finite homogeneous `H4` constrained
linear response; nonhomogeneous propagation, a physical tick, `c`, `G`,
Planck units and particle physics remain **NOT ESTABLISHED**.

