# Protocol: logical correction of the shifted augmented-rank classifier

Date: 2026-08-19

This is a fully target-disclosed classifier correction.  The direct audit's
literal `OPEN` outcome and all direct residual values are already known.  No
new physical calculation or blind claim is possible here.

## Frozen inputs

| input | SHA-256 |
|---|---|
| direct adversarial verifier | `72046f30f83e3af5192f60b108ea61e6b237dbecd11d615687b4c7c73417f521` |
| direct adversarial artifact | `9e9f7253fd10422f3534914fae020857162862123fd4eae889e3570083552179` |
| disclosed OPEN analysis | `c97c72260df77bb0509a35cd815d6a68147afb02e01f9f711d35af125fc068b3` |
| primary shifted verifier | `e4c5bcc18007c1c0ba7fbd38e29dffcc33a526fd790dbfcba8defe2ae44b7ab2` |
| primary shifted artifact | `0480f5d49d24e0f5d8e4e95f0cf62b7d0d9242459ed2b8f6d8e835ecd6e103a7` |
| direct adversarial protocol | `1dc9712a46b6ff6ac3c9b62e9d144f959f85622fe1ddbe2fc84de6ece3fa0982` |

The correction must preserve the direct artifact's
`SHIFTED_PSEUDOLONGITUDINAL_DIRECT_OPEN` outcome verbatim.

## Exact lemma

Let `X` and `Y` have the same row space carrier and let `P_X` be the
orthogonal projector onto `im X`.  Then

```text
rank([X,Y]) = rank(X) + rank((I-P_X)Y).
```

Proof: decompose every column of `Y` uniquely into its projection in `im X`
and its orthogonal component.  The projected part adds no dimension to
`im X`; the orthogonal components have trivial intersection with `im X`.

For the audit, `X=B L`, `Y=A L`, `rank(X)=15` and the span residual is
`||(I-P_X)Y||_2`.  Therefore a resolved nonzero span residual is already a
resolved witness that the exact augmented image is larger than the
15-dimensional image of `BL`.  A second numerical singular-value threshold
is a redundant sufficient test, not an independent mathematical hypothesis.

## Mechanical audit

The registered verifier will:

1. enforce every frozen hash;
2. require the direct artifact to remain `OPEN`, `18/18`, with 16 complete
   `15+10` stiffness cells;
3. inspect the frozen direct source and require that `left_bl` is constructed
   by full SVD of `BL`, while `span_residual` is the norm of
   `(I-left_bl left_bl*) AL`;
4. require in every direct cell `rank_longitudinal=15`, positive kinetic form,
   both denominators resolved, both residual inequalities true, and
   `rho_span` plus `rho_comm` labelled `NONZERO_RESOLVED`;
5. preserve and report the auxiliary failure: augmented classifier rank 15,
   with the sixteenth singular value positive but below its frozen threshold;
6. require the mechanically independent primary artifact to give the same
   two nonzero labels in all 16 cells;
7. classify only the persistence of non-invariance, not a numerical value for
   the redundant augmented rank.

## Outcome hierarchy

1. `SHIFTED_PSEUDOLONGITUDINAL_RANK_CORRECTION_CONTROL_FAILED` if any frozen
   provenance, source-shape, carrier or upstream control fails.
2. `SHIFTED_PSEUDOLONGITUDINAL_DIRECT_RESIDUAL_CONFIRMATION` if all 16 direct
   cells have resolved nonzero span and commutator residuals and both
   independent routes agree on those labels.
3. `SHIFTED_PSEUDOLONGITUDINAL_DIRECT_RESIDUAL_REFUTATION` if all 16 direct
   span and commutator residuals are `ZERO_CONSISTENT`.
4. `SHIFTED_PSEUDOLONGITUDINAL_RANK_CORRECTION_OPEN` otherwise.

## Interpretation firewall

Confirmation corrects an overly conjunctive numerical classifier.  It does
not convert the finite result into a symbolic theorem about the physical
background, and it does not establish conservation, curvature scaling,
refinement, continuum symmetry, instability or propagation.  The original
OPEN artifact is never overwritten.

Only this small verifier and static registry guards may be run.

