# First complete bounded-root result and deterministic-artifact defect

Date: 2026-08-20

Fast-evaluator implementation commit: `b38c451`.

The first complete 120-attempt execution passed `12/12` and wrote artifact
SHA-256

```text
17eeb5cb669e21547737b1e541334f581a42f0cd08bb5a9e2521900e59006515
```

Its frozen outcome was

```text
REFINED_H4_NO_FINITE_ROOT_FOUND_OTHER.
```

All 72 main-box and 48 boundary-ladder attempts ran. No main endpoint met the
`1e-7` refinement gate, no finite positive root was accepted, and zero of 12
classes met the strict zero-lapse boundary-pattern definition. The best
preconditioned main residual was `1.2091e-5`; the best ladder residual was
`2.6470e-5`.

This is a complete first execution, but not yet the accepted reproducible
result. The fast-evaluator correction protocol mistakenly required both
per-attempt wall-clock times and byte-identical deterministic JSON. Because
wall time is nondeterministic, those requirements cannot both hold.

The narrow correction is to remove elapsed seconds from the scientific JSON
while retaining timing only in the console/result narrative. No action value,
endpoint, optimizer status, residual, bound, seed, branch diagnostic,
classification or outcome may change. The corrected verifier must then run
twice with byte-identical artifacts before consolidation.
