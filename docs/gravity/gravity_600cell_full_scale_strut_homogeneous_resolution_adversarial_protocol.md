# Protocol: adversarial direct-matrix replication of the homogeneous line

Date: 2026-08-20  
Status: **preregistered post-result adversarial replication**

## Disclosed target and independence boundary

The primary result in commit `3ee5c55` is fully disclosed:
`HOMOGENEOUS_WEAK_POLE_LINE_UNIQUE`.  This verifier is allowed to compare with
that result, but it must not use its symbolic proof or its analytic generator to
construct the candidate line.

The decisive construction must instead:

1. rerun the frozen full-action geometry at P160;
2. select the unique one-dimensional constant-overlap representation sector;
3. rebuild the complete 65-by-10 carrier block directly from the edge lists and
   raw scale/strut response formulae;
4. rebuild the 65-by-5 weak canonical lift by a new Arb strong-block solve;
5. form the 65-by-10 D and 65-by-15 K intersection matrices;
6. extract their candidate null vectors by deterministic normal-equation solves,
   not SVD, eigensolvers, the stored P100 candidate or the primary analytic
   formula.

Using the first nine D columns and the first fourteen K columns as the dependent
blocks, set the final component to one and solve the overdetermined least-squares
normal equations.  This choice is frozen before execution.  No pivot or deleted
column may be selected after inspecting the result.

## Rank and residual requirements

At P160 in both parities require:

- every one of the ten D single-column-deleted Gram determinants excludes zero;
- every one of the fifteen K single-column-deleted Gram determinants excludes
  zero;
- the normalized direct-solve residuals
  `norm(M v)/(norm(M) norm(v))` are below `1e-30`;
- the five D scale components and five D strut components each have relative
  spread below `1e-30`;
- the K vector agrees, as a normalized projector, with the D vector joined by
  repeating its strut block, within `1e-30`;
- the even/odd D projectors agree within `1e-60`.

The nonzero deleted Gram determinants supply rank lower bounds 9 and 14.  The
accepted direct null residuals supply the complementary upper bounds in the
same frozen numerical model, hence one-dimensional kernels.  The exactness
claim still rests on the primary structural identity; this verifier attacks
whether that identity is the line present in the independently reconstructed
full matrices.

## Comparison and hostile controls

Only after extracting the direct D line, undo the already frozen carrier column
scales and compare its physical `sigma/c` with the primary value.  Require
absolute disagreement below `1e-30`.  This is a comparison, not a construction.

Controls:

- replacing `sigma=-lambda*p_z` by the missing-lambda direction must give a
  normalized D residual above `1e-10`;
- reversing the relative sign must give a normalized D residual above `1e-3`;
- a synthetic planted rank-nine 65-by-10 matrix must be accepted and an identity
  rank-ten control must be rejected as a one-line matrix;
- all pinned inputs must retain their frozen SHA-256 digests.

## Outcome hierarchy

1. `HOMOGENEOUS_ADVERSARIAL_CONTROL_FAILED`;
2. `HOMOGENEOUS_ADVERSARIAL_RANK_DISAGREEMENT`;
3. `HOMOGENEOUS_ADVERSARIAL_LINE_DISAGREEMENT`;
4. `HOMOGENEOUS_WEAK_POLE_LINE_REPLICATED`.

Outcome 4 upgrades the line to **DERIVED COMPUTATIONAL, adversarially
replicated**.  It still does not evaluate the omitted pole equation and does not
identify the line as gauge, dynamics, a clock or a physical tick.

Only this targeted verifier and the static registry guard may run.  No full
suite.

