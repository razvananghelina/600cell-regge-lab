# First canonical-data carrier failure

Date: 2026-08-19

## Complete hypotheses

This record concerns the fixed lower regular 600-cell, a flat Lorentzian
tetrahedral-frustum slab, variable face connections, and the frozen
target-blind augmented compatibility census in commit `e3a77fe`.  The tested
candidate carrier has one radial/scale variable and one normal/lapse variable
at each of the 120 vertices.  Its local displacement is

\[
\delta q_i=\sigma_v p_i+\nu_v n,
\]

and its induced canonical data are the derived upper-edge variation
`8 lambda (sigma_u + sigma_v)` and strut variation
`6 (lambda - 1) sigma_v - 2 tau nu_v`.  No action, Hessian, propagation speed,
or physical unit is evaluated here.

## Preserved first execution

The first registered execution produced
[`gravity_600cell_canonical_data_carrier.json`](../../reproducible/gravity_600cell_canonical_data_carrier.json)
with SHA-256
`e81d9b865a6a872fa894a8ccf68b718a38d122f31e42dc22f2232ed83c63bbfe`.
It reported 8/10 checks and the literal outcome
`CANONICAL_DATA_CARRIER_CONTROL_FAILED`.

The scientific candidate failed exact inclusion on all 3600 face rows for
both representatives `(lambda,tau)=(2,5)` and `(3,11)`.  Rebuilding the
coordinate graph with the alternate exact right inverse gave the same 3600
nonzero rows and the same first residual.  For `(2,5)`, row 2 is

\[
-\frac{25}{37}\sigma_0+\frac{25}{37}\sigma_{33}
+\frac{5}{37}\nu_0-\frac{5}{37}\nu_{33}.
\]

For `(3,11)`, row 2 is

\[
-\frac{242}{359}\sigma_0+\frac{242}{359}\sigma_{33}
+\frac{44}{359}\nu_0-\frac{44}{359}\nu_{33}.
\]

The candidate itself has exact data rank 240.  The endpoint-difference and
deleted-lapse negative controls also fail inclusion, as required, and the two
constant columns reproduce the frozen homogeneous scale/lapse controls.

## Diagnosis

**DERIVED NEGATIVE.** The proposed
`Q^120 vertex scale + Q^120 vertex lapse` carrier is not the admissible
canonical-data kernel.  A dimension match `240 = 120 + 120` did not establish
inclusion, and the exact inclusion test refutes it.

**STRUCTURAL SOFTWARE DEFECT.** The literal `CONTROL_FAILED` label is caused by
the classifier treating successful alternate-graph *inclusion* as a control.
The actual construction control should require that the baseline and alternate
graphs agree on the inclusion decision.  Here they agree that inclusion is
false.  Correcting that decision hierarchy must not alter the matrices,
candidate map, exact residual criterion, or preserved first artifact.

**OPEN.** The frozen modular nullity 240 remains only a stable modular result
and an upper bound on the rational admissible dimension.  Its actual rational
carrier is not classified by this failure.
