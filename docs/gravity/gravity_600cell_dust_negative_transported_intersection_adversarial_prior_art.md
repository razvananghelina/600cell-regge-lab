# Independence gate: adversarial negative-intersection replication

Date: 2026-08-18

Status: **PRIMARY RESULT DISCLOSED; NO FLOAT64 AUDIT SPECTRUM INSPECTED.**

## Purpose

The preregistered high-precision calculation has reported

```text
NEGATIVE_TRANSPORTED_INTERSECTION_ZERO_CERTIFIED_ALL
```

with exact rank `30` in all `16` cells.  Two byte-identical executions prove
reproducibility, not independence.  This post-result gate attempts to falsify
the result using the earlier binary negative projectors, the earlier tangent
archive and a different square leakage matrix.

For orthonormal phase bases `W_0,W_1` and an orthonormal target complement
`W_1_perp`,

```text
L = W_1_perp* T_2 W_0 : C^30 -> C^30
```

has kernel equal to the transported intersection.  Hence full rank of `L` is
equivalent to zero intersection.

## Literature and scope

The SVD/principal-angle description and perturbation of subspace
intersections are standard; see:

- Knyazev and Argentati, *Principal Angles between Subspaces in an A-Based
  Scalar Product*, SIAM J. Sci. Comput. **23** (2002), DOI
  `10.1137/S1064827500377332`;
- Knyazev, Jujunashvili and Argentati, *Angles between infinite dimensional
  subspaces with applications to the Rayleigh-Ritz and alternating projectors
  methods*, J. Funct. Anal. **259** (2010), DOI
  `10.1016/j.jfa.2010.05.018`.

The physical restriction remains Dittrich and Hoehn's action-derived
pre/post-constraint formalism, DOI `10.1063/1.4818895`, arXiv:`1303.4294`.
The present Hilbert-metric negative fiber is not established to be such a
constraint surface.

None of these sources determines the `600`-cell rank.  External novelty is
**OPEN**.

## Mechanically different path

The primary calculation used:

- Flint source and tangent balls;
- `mpmath` at `100` digits;
- newly reconstructed high-precision projectors;
- the `60 x 60` residual `(I-Q_1)T_2Q_0` and its structural rank bound.

The audit will instead use:

- the binary projectors committed before the primary result;
- a fresh `scipy.linalg.eigh` extraction of their rank-`15` ranges;
- the earlier committed binary tangent midpoint archive;
- a complete-QR target complement;
- the square `30 x 30` leakage and binary LAPACK SVD.

It will also compare a separately constructed null-space complement and a
second SVD driver, reverse/rephase bases, and swap the `(q,p)` block-ordering
convention together with the tangent and fibers.

This route shares the frozen mathematical object but not the decisive
high-precision projector, tangent ball, residual matrix or SVD
implementation.

## Controls and limitation

For every cell:

- `T_full=W_1 W_0*` must give zero leakage and a `30`-dimensional
  intersection;
- `T_zero=W_1_perp W_0*` must give full leakage rank and zero intersection;
- basis gauge, QR/null-space complements, two SVD drivers and `(q,p)` ordering
  must agree within a preregistered roundoff floor.

The audit is binary64 and does not propagate the old projectors' physical
source enclosures.  It is therefore not a second exact certificate.  Its role
is **STRUCTURAL INDEPENDENT CORROBORATION** of a robust rank.  Any actual cell
which is not full rank, or any failed stress, returns the consolidated verdict
to **OPEN** under project rule 4.
