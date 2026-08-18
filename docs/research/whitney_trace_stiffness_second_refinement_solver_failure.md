# Recorded solver failure: second trace-stiffness refinement

Date: 2026-08-11

Protocol commit: `702fa5b`

## Status

The preregistered three-level verifier returned **10/11**.  Levels zero and
one reproduced all dense spectral edges to (7.11\times10^{-15}) relatively.
All exact topology, kernel, face-metric, and row-image checks passed.

At level two, the generalized LOBPCG blocks failed the unchanged
(10^{-7}) Ritz gate:

- degree zero maximum block residual: (4.32\times10^{-4});
- degree one maximum recorded residual: (8.39\times10^{-4}).

Degree two passed with residual (3.53\times10^{-12}).

## Provisional values are rejected

The run printed provisional step-`1->2` ratios

```text
(7.2761, 5.4103, 5.4488)
```

and provisional degree spread `1.34486`.  They are **NOT ACCEPTED** because
the complete preregistered extremal gate failed.  This commit preserves them
only as solver provenance.

## Structurally stronger method available

For the 10,980- and 13,860-dimensional level-two quotients, the two reduced
matrices can be formed explicitly and sparsely:

\[
 A=V^*HRM^{-1}R^*HV,
 \qquad
 B=V^*HV.
\]

Both are real symmetric and (B>0).  Direct sparse generalized symmetric
Lanczos can target their edges without changing the quotient, operator, or
residual gate.

That is a post-failure solver correction.  It must be frozen in a separate
commit and calibrated on levels zero and one before the provisional values
can be reconsidered.

## Ledger

- **DERIVED:** complete second-refined carrier and exact metric-type geometry.
- **DERIVED:** exact kernels and positive degree-two quotient.
- **FAILED GATE:** level-two generalized LOBPCG for degrees zero and one.
- **NOT A RESULT:** provisional step ratios and balance trend.
- **OPEN:** explicit sparse-quotient correction.
