# Correction protocol: generic-velocity leading verifier

Date: 2026-08-21

Registered implementation commit: `dfa2688`.

Preserved first-run artifact:

```text
reproducible/gravity_600cell_generic_velocity_composition_first_failure.json
SHA-256 feab4b4a66422df04ecc5906fcb6bbec016365ea1766e9f6810afa8aec49a275
```

The first run returned `8/11` and outcome
`GENERIC_VELOCITY_LEADING_OPEN`.  All action, lapse, momentum,
time-reversal, composition and 100-decimal controls passed.  The failures are
resolved as follows before any rerun.

## 1. Rational-form equality

The positivity-domain gate compared the structurally different but equal
forms

```text
(v^2+2)/(2 v^2+6)
and
(v^2+2)/(2(v^2+3))
```

using Python structural equality after `factor`.  Replace both bound checks
and the substituted constraint residual by

```text
simplify(left-right)==0.
```

No domain, branch or inequality changes.

## 2. Correct static mass curvature

The verifier independently derived

```text
theta'(v)=-2v/[(v^2+3)sqrt(v^2+4)sqrt(3v^2+8)],
```

so

```text
theta''(0)=-1/(6 sqrt(2)),
epsilon_v''(0)=5/(6 sqrt(2)).
```

The hand-entered comparison target mistakenly used twice this value.  Freeze
the corrected consequence

```text
mu''(0)
 = (180/pi)[5/(12 sqrt(2))-epsilon/8]
 = 15[5 sqrt(2)-3 epsilon]/(2 pi)
 > 0.
```

This coefficient was not a preregistered target; the protocol required the
verifier to derive it and determine its sign.  The computed derivative is
retained and the incorrect hand check is corrected to match it.  Also record
the individual branch/static booleans in the next artifact.

Do not change the action, branch count, velocity/mass domain, numerical
points, precision, convergence window or outcome hierarchy.  Rerun only the
targeted verifier.

