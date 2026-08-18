# Result: the canonical kinetic norm does not close the safety margin

Date: 2026-08-18

## Headline

**DERIVED COMPUTATIONAL, TARGET-DISCLOSED.**  The action-selected kinetic
metric is positive-resolved on every one of the `16` frozen negative-shape
cells.  Its unique positive square root defines the preregistered canonical
similarity, and every literal matrix-valued Rouche cover still passes.

All `16` covers fail the repository's separate `100 x error` rule.  The
preregistered outcome is therefore

```text
NEGATIVE_SHAPE_KINETIC_SAFETY_OPEN.
```

The canonical norm does not strengthen the previous Euclidean certificate.
Its smallest sampled signal/error ratio is approximately `4.96327`, compared
with `8.16755` in the Euclidean norm.  This closes the numerical-norm branch:
trying another similarity after seeing these results would be fitting.

## Provenance and reproduction

| stage | commit |
|---|---|
| primary-literature and structural gate | `cdcaf8d` |
| target-disclosed protocol | `afe2c4e` |
| frozen carrier export | `6bd1a02` |
| verifier registered before first execution | `c232291` |

Verifier:

```text
reproducible/verify_gravity_600cell_dust_negative_shape_kinetic_metric.py
```

Artifact:

```text
reproducible/gravity_600cell_dust_negative_shape_kinetic_metric.json
SHA-256 22a6ab9e45179773f84505b1ee6fd9bf940a70c6c9d08cfd1b55c85d490c6ed7
```

Two targeted executions were byte-identical and each reported `10/10`.  The
full suite was deliberately not run, following the user's instruction.

## What was tested

For every inherited negative carrier, the verifier uses

```text
B_- = E_-^* (-M_S) E_-,
S   = B_-^(1/2),
Q_B(z) = S Q(z) S^(-1).
```

The source coefficient balls are transported conservatively:

```text
epsilon_Gamma,B = kappa_2(S) epsilon_Gamma,
epsilon_Omega,B = kappa_2(S) epsilon_Omega.
```

No eigenvector balancing, diagonal optimization or selected angular weighting
enters the calculation.  The exact roots are unchanged by similarity; only
the sufficient singular-value bound is being audited in the inherited norm.

## Numerical ledger

Across the `16` cells:

```text
lambda_min(B_-) / epsilon_B   352.152279211 ... 352.152279223
kappa_2(B_-)                   2.70824810972833 ... 2.70824810972834
kappa_2(S)                     1.64567557851732 ... 1.64567557851733
literal weakest sampled ratio  4.96326612936 ... 4.96326613280
literal evaluated intervals   18832 in every cell
literal maximum depth             16 in every cell
```

The weakest certified literal leaf lower bound is between
`6.7791501e-10` and `6.7791572e-10`.  The first `100x` failure occurs at
`theta = 0.00306796157577`; its unscaled signal/error ratio is approximately
`45.3525`, already below `100`.  As in the Euclidean audit, that is the first
failure in traversal order, not the global minimum.  The completed literal
cover later samples the smaller ratio `4.96327`.

The square-root controls are comfortably inside their numerical floors.  The
largest ratios to the corresponding control floors were approximately:

```text
S^2 reconstruction                 1.19e-3
S S^(-1) reconstruction            5.92e-6
scalar-metric rescaling matrix      2.04e-3
scalar-metric rescaling condition   1.63e-4
unitary reversal singular values    1.84e-3
```

Multiplying `B_-` by the preregistered scalar `7` changes no cover verdict,
and the fixed unitary reversal preserves the sampled singular values.

## Hostile interpretation audit

1. **This is a negative about a sufficient bound, not about the exact roots.**
   The inherited literal coefficient balls still have the previously
   transferred root count `15/0/15`; the stronger project margin remains
   unresolved.
2. **The kinetic norm is canonical only relative to the frozen local action.**
   It is selected independently of the root target, but it does not prove that
   independently solved later slabs share the same recurrence.
3. **A different norm is no longer an honest continuation.**  The Euclidean
   and the only action-selected kinetic norm have now both been evaluated.
   Choosing another similarity to improve the observed margin would introduce
   post-result fitting.
4. **No physical instability is yet established.**  A frozen `15/0/15` split
   is not a Lyapunov exponent, a polarization identification, a propagation
   speed or nonlinear long-time growth.
5. **The result does not damage the positive kinetic form.**  Positivity is
   resolved in all cells and is not the source of the failed `100x` gate.

## Status ledger

- **DERIVED COMPUTATIONAL:** all `16` inherited kinetic restrictions are
  positive-resolved by more than `352` error units.
- **DERIVED COMPUTATIONAL:** all `16` literal kinetic-norm Rouche covers pass.
- **DERIVED UPSTREAM:** their transferred exact root count is `15/0/15`.
- **OPEN BY PREREGISTERED SAFETY RULE:** all `16` kinetic-norm `100x` covers
  fail; the worst sampled literal ratio is only `4.96327`.
- **REFUTED AS A ROUTE:** the action-selected kinetic norm does not upgrade the
  frozen root split to the project's `100x` physical-resolution convention.
- **OPEN:** nonlinear carrier invariance, non-autonomous persistence,
  polarization, continuum dispersion, refinement and physical growth.

## Next load-bearing gate

Stop optimizing the frozen polynomial.  The next admissible physical mission
is to solve two or more nonidentical consecutive slabs independently, derive
their time-dependent negative-shape recurrence from those solutions, and ask
whether the same carrier and root separation persist.  Its protocol must fix
the slab family, continuation variables and matching conditions before any
growth spectrum is inspected.

