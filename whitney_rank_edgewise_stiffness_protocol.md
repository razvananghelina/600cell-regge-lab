# Protocol: Whitney trace stiffness on the rank-edgewise tower

Date: 2026-08-11

This protocol is frozen before constructing or diagonalizing the new
carriers.  No phenomenological target, fitted exponent, preferred constant,
or favorable subset of form degrees will be used.

## Fixed question

Does the degree balancing previously seen under iterated barycentric
refinement survive when the carrier is replaced by the canonical,
uniformly shape-regular tower

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}K),
 \qquad k\in\{1,2,4\}?
\]

The closed control `K` is the boundary of the 4-simplex, with five regular
tetrahedra.  This is the same calibrated control used in the existing exact
trace-stiffness calculations.  The operator is unchanged: exact local Whitney
mass, exact shared-face trace mass, the full jump map, and the quotient by its
exact row image.

## Carrier construction

1. Apply the already-certified barycentric subdivision once.
2. Its chamber vertices are ordered by face rank `(0,1,2,3)`.
3. Construct `Esd_k` directly from exact color schemes for `k=1,2,4`.
4. Merge a new vertex globally by its exact rational barycentric weights on
   barycentric vertices; do not merge by numerical tolerance.
5. Retain every cell and every exact local metric type.

The expected top-cell counts are fixed in advance:

\[
 120k^3=(120,960,7680).
\]

Lower-dimensional f-vectors, metric-type counts, gaps, and ratios are outputs,
not targets.

## Numerical and exact gates

- Every shared face must have exactly the same induced metric from both
  parents.
- Every occurrence graph must be connected and every quotient rank must equal
  the exact combinatorial rank.
- The existing `k=1` barycentric certificate must be reproduced within
  relative error `5e-7` for every positive gap and maximum eigenvalue.
- Every reported Ritz residual must be below `1e-7`.
- Because the normalized shape set is identical for `k=2` and `k=4`, the
  worst local Dirac norm must scale by exactly two, checked numerically to
  relative error `1e-10`.

If the `k=4` sparse solve exceeds available resources, that is a recorded
solver boundary, not a numerical result.  The `k=1 -> 2` step alone cannot
establish or refute repeated scaling.

## Frozen comparisons

For form degree `p=0,1,2`, define

\[
 s_{k,p}=a_k/g_{k,p},
 \qquad
 R^{12}_p=s_{2,p}/s_{1,p},
 \qquad
 R^{24}_p=s_{4,p}/s_{2,p},
\]

where `a_k` is the worst local Kähler--Dirac norm and `g_{k,p}` is the
smallest positive trace-stiffness quotient eigenvalue.

Report all six ratios.  Define the degree spread of a ratio vector by

\[
 \operatorname{spread}(R)=\max_pR_p/\min_pR_p.
\]

- **DERIVED NUMERICAL exact repetition** only if every component of
  `R24/R12` differs from one by less than `1e-6`.
- **PATTERN toward degree balance** only if
  `spread(R24) < spread(R12)`.
- **PATTERN NEGATIVE** if that spread does not improve.
- Equality of two selected degrees is never tested separately and carries no
  label.

The comparison will also be made with the two-step barycentric control, but
only after the new ratios are fixed.  That comparison is diagnostic and
cannot change the labels above.

## Acceptance and kill boundaries

**Acceptance for the scaling clue:** the exact construction and solver gates
pass, and the second-step degree spread is smaller than the first-step spread.
This would be a pattern worth testing on the complete 600-cell, not a physical
renormalization law.

**Kill for the clue on this control:** exact construction passes but degree
spread is unchanged or worse.  Then the earlier apparent balancing was caused
by, or depended on, barycentric shape degeneration.

Regardless of outcome, this experiment cannot select an absolute stiffness,
length, tick, mass, causal speed, or Planck scale.

