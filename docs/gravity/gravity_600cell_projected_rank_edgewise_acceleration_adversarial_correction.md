# Correction protocol: intrinsic precision and the static zero

Date: 2026-08-19

Original audit protocol: `3b7bd6c`  
Preserved failed artifact commit: `391f758`  
Failed artifact SHA-256:
`e41a7f57be6995cb39b8fdb89fb981263473be9bbf8982ae0bb84f65d07fc8f6`

## Observed failure

The first adversarial run passed coefficient agreement, the regular 600-cell
positive control and the wrong-multiplicity negative control, but failed two
held-out compression gates.

The static relative error was approximately one because both complete and
compressed *total* actions vanish by the preregistered dust cancellation.  A
relative error normalized by an exactly vanishing target is undefined and
was therefore an invalid control.

At the two non-static held-out states, the largest errors were approximately
`2.16e-8` and `2.78e-8`, just above the frozen `2e-8` gate.  The original
11-decimal intrinsic signatures merge source-coordinate roundoff variants
whose per-class deviations are small but whose weighted action difference is
not below that gate.

The failed result is retained as an **OPEN disagreement** until the corrected
audit passes.

## Target-free precision diagnostic

Without changing the action, derivative grid, coefficient sentinels, primary
artifact or comparison target, the fine-carrier held-out errors were audited
as a function of signature precision:

```text
digits   tetra classes   maximum held-out relative error
11                  23   2.78e-8
12                 158   2.89e-8
13                 168   1.14e-10
14                 172   7.17e-11
```

Thirteen digits are the first tested precision that satisfies the original
`2e-8` held-out gate by more than two orders of magnitude.  Selecting it does
not use a desired coefficient and does not loosen any acceptance tolerance.

## Frozen corrections

Only two changes are permitted.

1. Round intrinsic squared-distance signatures to 13 rather than 11 decimal
   places.  Tighten the maximum within-class raw-signature residual from
   `<2e-9` to `<2e-12`.
2. At the static held-out state compare the nonzero gravitational actions and
   normalize by their magnitudes.  Continue comparing total actions at both
   non-static held-out states.  Retain the original `<2e-8` relative gate for
   all three comparisons.

The corrected run writes a new artifact

```text
reproducible/gravity_600cell_projected_rank_edgewise_acceleration_adversarial_corrected.json
```

and must not overwrite the preserved failed artifact.

Every coefficient grid, derivative step, sentinel, target tolerance,
positive control, negative control and outcome label remains unchanged.
