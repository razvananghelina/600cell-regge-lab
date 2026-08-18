# Framing gate: resolve the adversarial negative-intersection roundoff

Date: 2026-08-18

Status: **PRIMARY AND FIRST AUDIT RESULTS DISCLOSED; NO FIXED-INPUT
HIGH-PRECISION LEAKAGE SPECTRUM COMPUTED.**

## The disagreement that must be resolved

The source-certified primary verifier reported exact leakage rank `30` and
zero transported intersection in all `16` cells.  The mechanically different
binary64 audit reproduced the same qualitative singular spectrum but its
preregistered roundoff floor

```text
1000 eps_machine 60 max(1,||T_2||,||L||)
```

was larger than the final directions.  It therefore resolved rank only `20`
in every cell and returned

```text
ADVERSARIAL_NEGATIVE_INTERSECTION_DISAGREEMENT_OPEN.
```

That outcome is retained.  Its threshold may not be changed retroactively.
Under project rule 4 the consolidated scientific result is presently
**OPEN**.

## Exact new question

Treat every stored binary64 entry of the earlier negative projectors and
tangent midpoint as an exact dyadic input.  Does the corresponding square
leakage matrix

```text
L_bin = W_1,perp* T_bin W_0
```

have rank `30` when all eigenspaces, products and singular values are
recomputed at independently varied high precision with residual-based error
bounds?

This question distinguishes:

1. a genuine rank loss already present in the independent binary data; from
2. a conservative binary64 arithmetic floor which was too broad to classify
   small but genuine singular values.

It does **not** replace the primary Flint source enclosure.  Rank of an exact
stored midpoint does not prove rank of the underlying source ball.

## Numerical basis and literature

For a Hermitian invariant-subspace residual separated by a spectral gap,
Davis--Kahan-type bounds control the projector rotation.  Singular values then
obey Weyl's perturbation inequality.  The relevant primary numerical sources
remain:

- Knyazev and Argentati, DOI `10.1137/S1064827500377332`;
- Knyazev, Jujunashvili and Argentati, DOI
  `10.1016/j.jfa.2010.05.018`.

The physical pre/post-constraint caveat remains Dittrich and Hoehn, DOI
`10.1063/1.4818895`, arXiv:`1303.4294`.

These sources supply no `600`-cell rank.  External novelty remains **OPEN**.

## Mechanically distinct fixed-input certificate

The calculation will not reuse the primary Flint `M,V` balls or its negative
projectors.  It starts from the earlier binary projectors and tangent archive,
converts each stored float to its exact dyadic rational through
`float.as_integer_ratio()`, and then:

- diagonalizes each binary projector at `100` and `140` decimal digits;
- uses top/bottom spectral ranges directly as phase fiber/complement;
- bounds their invariant-subspace error from cross residual and spectral gap;
- forms the square `30 x 30` leakage;
- computes singular values both by direct SVD and through the Gram spectrum;
- checks the smallest value independently through `1/||L^-1||`;
- requires complete two-precision agreement and known rank-zero/rank-30
  controls.

The primary used source balls, stiffness eigenspaces and a `60 x 60` residual.
The first audit used SciPy/NumPy, complete QR and a global float64 floor.  This
bridge uses fixed binary inputs, high-precision projector eigenspaces, their
spectral complements and a square leakage.  It is a resolution of the
arithmetic disagreement, not a third discovery attempt.

## Acceptance logic

- If every fixed-input cell has all `30` singular values separated from the
  complete residual error at both precisions, the first audit's open label is
  attributed specifically to its preregistered global roundoff envelope.  In
  combination with the independent binary spectrum and the separate Flint
  source certificate, project rule 4 is satisfied.
- If even one cell fails, remains precision-sensitive or violates an
  SVD/Gram/inverse/control check, the consolidated result remains **OPEN**.

No empirical spread from the disclosed audit is used as an error bar, and no
observed singular value is used to tune a threshold.
