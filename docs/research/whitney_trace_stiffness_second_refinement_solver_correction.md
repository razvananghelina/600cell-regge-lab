# Post-failure correction: explicit sparse trace quotient

Date: 2026-08-11

Original protocol: `702fa5b`  
Recorded failure: `6f21619`

## Correction frozen before recomputation

For degrees zero and one, form the reduced quotient matrices explicitly:

\[
 A=V^*HRM^{-1}R^*HV,
 \qquad
 B=V^*HV.
\]

The factors are the same exact trace jump (R), face metric (H), local
element mass (M), and row-image basis (V) used by the failed matrix-free
LOBPCG.  No operator or weight changes.

The matrices must be sparse and numerically symmetric, with maximum transpose
residual below (10^{-11}).  Their dimensions must equal the exact quotient
ranks.

Use generalized symmetric Lanczos (`eigsh`) with:

- five smallest eigenpairs through shift-invert `sigma=0`, `which="LM"`;
- five largest eigenpairs with `which="LA"`;
- tolerance (10^{-11});
- maximum 20,000 iterations;
- the original deterministic seed rule;
- unchanged directly recomputed generalized Ritz gate (10^{-7}).

Degree two keeps the already corrected full-row symmetric Lanczos method.

## Mandatory calibration and cross-checks

1. The explicit method must reproduce all degree-zero and degree-one dense
   edges at levels zero and one within the original
   (5\times10^{-7}) relative gate.
2. It must reproduce the accepted matrix-free level-one quotient edges within
   the same gate.
3. Every level-two smallest and largest block must pass the unchanged
   (10^{-7}) Ritz gate.

Only after all three pass may the step-`1->2` ratios be interpreted.  The
provisional values from `6f21619` are not targets and are not acceptance
criteria.
