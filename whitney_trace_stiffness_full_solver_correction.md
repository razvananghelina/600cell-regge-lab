# Post-failure solver correction: full-row degree-two trace quotient

Date: 2026-08-11

Original protocol: `7ed6d49`  
Recorded failure: `9c04a94`

## Frozen correction before recomputation

The original generalized block LOBPCG remains unchanged for degrees zero and
one.  Their complete recomputed Ritz residuals already pass the frozen
(10^{-7}) gate.

Only degree two uses the following exact structural reduction.  Each shared
triangle contributes one 2-form jump, every triangle has two tetrahedral
parents, and the degree-two occurrence graph has two nodes and one edge.
Consequently:

\[
 \operatorname{rank}R_2=\#\text{rows}(R_2),
\]

so there is no row-cycle nullspace and no generalized row-image metric to
solve.

The positive spectrum is exactly that of

\[
 K_2=H_2^{1/2}R_2M_2^{-1}R_2^*H_2^{1/2},
\]

an ordinary real symmetric positive sparse operator.

Use symmetric Lanczos (`scipy.sparse.linalg.eigsh`) with:

- five smallest eigenpairs using `which="SA"`;
- five largest eigenpairs using `which="LA"`;
- tolerance (10^{-11});
- maximum 20,000 iterations;
- deterministic initial vector from the same frozen seed rule;
- unchanged recomputed relative Ritz gate (10^{-7}).

The square root (H_2^{1/2}) is unambiguous because every degree-two face
mass block is a positive scalar.

## Mandatory checks

1. Reproduce both dense degree-two control gaps and maxima within the original
   (5\times10^{-7}) relative calibration gate.
2. Reproduce the complete base degree-two LOBPCG edge within that gate.
3. Certify the refined smallest and largest five blocks below the unchanged
   (10^{-7}) Ritz gate.

No physical target, provisional ratio, or desired spread enters this repair.
If any check fails, the complete trace result remains OPEN.
