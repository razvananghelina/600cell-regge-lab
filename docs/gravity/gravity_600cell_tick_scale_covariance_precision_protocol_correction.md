# Pre-evaluation correction: remove an invalid cross-parity action equality

Date: 2026-08-21

Original precision protocol: commit `093055e`.

Status: frozen before implementation and before any arbitrary-precision direct
2400-simplex evaluation.

## Error

Section 4 of the original protocol required the even and odd high-precision
base actions to agree within `1e-55`.  That is not a legitimate control for
the frozen off-shell state.  Its deterministic perturbation is applied by
orbit-coordinate index separately in the two parity bases.  Those bases are
not identified by a frozen physical edge transport, so the two perturbed
vectors need not describe the same labelled metric.

This was detectable from the already committed primary artifact, whose base
actions are

```text
even: 3.4055269496999591195...e-6,
odd : 3.5879031191437823507...e-6.
```

No precision-adjudication value was evaluated before this correction.

## Corrected gate

Delete only the cross-parity base-action equality.  Retain every within-parity
scale identity and require each direct high-precision base action to agree
with its *own* stored primary orbit action within `1e-45`.

No parity averaging, relabelling, fitted transport, tolerance change or action
change is admitted.  All other hypotheses, evaluations, gates and outcomes of
commit `093055e` remain frozen.

