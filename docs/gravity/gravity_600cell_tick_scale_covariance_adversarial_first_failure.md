# First adversarial result: covariance passes, cross-implementation gate fails

Date: 2026-08-21

## Provenance

- frozen protocol: `fd1f8a8`;
- primary result: `0ac0aba`;
- registered adversarial implementation before evaluation: `b94d7ec`;
- targeted verifier:
  `reproducible/verify_gravity_600cell_tick_scale_covariance_adversarial.py`;
- artifact:
  `reproducible/gravity_600cell_tick_scale_covariance_adversarial.json`;
- artifact SHA-256:
  `9a28debdade3f189e3c21416af3a21a4516f7aa49f0686b3b39952ac87774d2e`.

Only the new targeted verifier was run.  It returned **10/12** and exited
nonzero.  The failure is retained rather than hidden by changing a threshold.

## Mechanical outcome

```text
TICK_SCALE_COVARIANCE_IMPLEMENTATIONS_DISAGREE
```

## What passed

The mechanically different binary64 implementation sums all 2400
four-simplices directly and differentiates with respect to raw squared
lengths.  For both parities and both frozen scale factors:

```text
S_scaled = alpha^2 S_base,
dS_scaled/dq = dS_base/dq
```

within the preregistered propagated bounds.  The largest observed action,
current-boundary-gradient and old-boundary-gradient errors are respectively

```text
6.575e-8,
1.350e-7,
1.154e-10,
```

while the propagated bounds lie between `1.311e-6` and `3.642e-6`.
All 35 internal and all 60 boundary raw-gradient components are nonzero at the
frozen support threshold.  All fixed-mass hostile controls fail covariance by
resolved amounts.

## What failed

The direct binary64 base action disagrees with the 100-decimal orbit action by

```text
even: 1.648e-6,
odd : 1.338e-6.
```

The implementation imposed a separate inherited comparison gate of `2e-8`.
That gate is tighter than the direct evaluator's own conditioning-aware
binary64 envelopes.  Because it was frozen before evaluation, it cannot be
relaxed after seeing these numbers.

This is not evidence against scale covariance: each implementation separately
obeys the scaling identity, including all derivative components.  It is also
not acceptable evidence for agreement: the required absolute cross-check did
not pass.

## Honest status

- **DERIVED COMPUTATIONAL:** direct covariance and its negative control pass.
- **OPEN:** whether the base-action discrepancy is binary64 loss near the thin
  Lorentzian slab or a genuine orbit/direct construction mismatch.
- **NOT ACCEPTED:** the absolute classical tick no-go remains provisional
  despite the exact dimensional proof.

## Frozen repair direction

Do not change the failed verifier.  Preregister a separate adjudication that
reconstructs the action by a direct 2400-simplex sum at arbitrary precision,
without orbit multiplicities, and compares it independently to both stored
values.  If high precision selects neither implementation, the disagreement is
the result.  If it selects the primary orbit action and also obeys the scale
identity, the binary64 absolute discrepancy is a resolved precision limit.

