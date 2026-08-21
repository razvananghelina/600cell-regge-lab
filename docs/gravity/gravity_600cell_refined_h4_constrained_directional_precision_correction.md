# Preregistered correction: high-precision direction reconstruction

Date: 2026-08-21

The default-precision control failure is preserved in commit `136da68`.

Move only the following frozen-artifact conversions and derived matrix
products into one `mp.workdps(180)` context:

```text
primary boundary basis,
24 stored primary lifts and response matrices used by the 12 directions,
stored quadratic values,
full boundary+internal directions,
coarse/finest displacement values.
```

The resulting 180-digit objects may then be rounded by each 140/180-digit
action-evaluation context.  Also perform the post-ladder error, ratio,
envelope, polynomial and corruption arithmetic inside `mp.workdps(180)`;
otherwise a second default-precision truncation would remain.  Require the
existing first-Richardson reproduction error `<1e-45`; it must falsify an
incomplete repair.

Do not change the stored strings, action, steps, extrapolation, ratio interval,
envelopes, polynomial or corruption controls, outcome hierarchy or scope.
Rerun only this diagnostic twice and require a byte-identical artifact.
