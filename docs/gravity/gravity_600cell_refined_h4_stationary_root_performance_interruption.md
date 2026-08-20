# Performance interruption: the literal 50-decimal `3-point` search is impractical

Date: 2026-08-20

Lorentzian-gate implementation commit: `c8d21b0`.

After all eight pre-search controls passed, the verifier began the 72 main-box
attempts. The process consumed one CPU core continuously for approximately 20
minutes without finishing the first six-seed class and was interrupted with
`SIGINT` inside SciPy's numerical-Jacobian evaluation.

No main or ladder artifact was written and no endpoint or scientific outcome
was classified. The previous `CONTROL_FAILED` JSON remains the last complete
artifact.

The bottleneck is mechanical: `jac='3-point'` evaluates the 50-decimal complex
Regge action roughly twenty additional times per optimizer Jacobian, while
SciPy's `nfev` does not expose that full cost in the progress log. Extrapolating
this incomplete class would put the 120-attempt census in the multi-hour to
multi-day range without checkpointing.

A repair may replace only the search-phase evaluator by a mechanically
equivalent complex-binary64 implementation if it is preregistered and compared
componentwise against the frozen 50/80-decimal evaluator at branch-valid
anchors. The 120 attempts, `3-point` solver, bounds, seeds, residual
preconditioner and high-precision candidate validation must remain unchanged.
Failure of the cross-evaluator control must stop the search.
