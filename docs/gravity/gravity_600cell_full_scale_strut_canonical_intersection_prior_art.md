# Prior-art gate: complete scale--strut/canonical intersection

Date: 2026-08-20

## Exact question and hypotheses

Fix the accepted non-static Lorentzian 600-cell Regge--dust slab, separately
for both frozen staircase parities.  Require all of the following:

```text
the regular 600-cell lower and upper boundary carriers;
the accepted curved background and its frozen dust mass;
real lambda != 1 and tau != 0;
(lambda-1)^2-3 tau^2 != 0;
the exact rank-240 geometry-selected scale--strut carrier G;
the same 840 internal plus 720 new-boundary log-squared-length row order;
the frozen Regge-action Hessian and its canonical strong-equation graph C;
the free binary-tetrahedral 2T action and its seven minimal sectors;
no extra face connection, Schur coefficient or continuum mode target.
```

In a minimal sector of irrep dimension `d`, write

```text
G = [G_scale,G_strut] : C^(10d) -> C^(65d),
C                       : C^(5d)  -> C^(65d).
```

Both maps use the same five pole-orbit positions.  Geometry gives the pole
blocks

```text
G_scale|pole = 0,
G_strut|pole = I,
C|pole       = I.
```

Therefore a common image vector must have canonical coefficients equal to
the strut coefficients.  The complete intersection is exactly parameterized
by

```text
ker D,  D = [G_scale, G_strut-C] : C^(10d) -> C^(65d).
```

This reduction is valid only after all three literal pole identities and
full-column-rank inputs are checked.  Independently,

```text
dim(im G intersect im C) = rank(G)+rank(C)-rank([G,C]).
```

must give the same number.

## Primary-literature boundary

- Barrett et al., [*A Parallelizable Implicit Evolution Scheme for Regge
  Calculus*](https://arxiv.org/abs/gr-qc/9411008), evolve vertices locally and
  apply the scheme to a 600-cell Friedmann cosmology.
- De Felice and Fabri study dust evolution and generalized 600-cell
  evolution in [arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093)
  and [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077).
- Dittrich and Hoehn, [*From covariant to canonical formulations of discrete
  gravity*](https://arxiv.org/abs/0912.1817) and [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974), derive action-generated
  pre/post constraints and canonical discrete evolution.
- Hoehn, [*Canonical linearized Regge Calculus: counting lattice gravitons
  with Pachner moves*](https://arxiv.org/abs/1411.5672), identifies
  curvature-carrying modes only after Hessian and gauge analysis on flat
  backgrounds.
- Tsuda and Fujiwara, [*Oscillating 4-Polytopal Universe in Regge
  Calculus*](https://arxiv.org/abs/2011.04120), derive symmetry-reduced
  4-polytopal FLRW dynamics and refinement limits.

These works make 600-cell evolution, action-generated canonical maps and
Hessian mode analysis **KNOWN**.  None of the located primary sources gives
this exact nonhomogeneous `120+120` carrier or its intersection with this
dust-slab strong-equation graph.  Search absence is not a novelty proof;
external novelty remains **OPEN**.

## KNOWN / CONTROL / OPEN

- **DERIVED INPUT:** `G` is geometry-selected and has exact global rank 240.
- **DERIVED INPUT:** `C` is action-generated and has rank 120 in all frozen
  sectors.
- **DERIVED NEGATIVE INPUT:** the pure-strut half alone has zero intersection
  with `C` in all 14 parity/sector cases.
- **CONTROL:** the pole identities make the reduced-kernel formula exact.
- **CONTROL:** the unreduced joined-image formula must agree.
- **OPEN:** every actual full-carrier intersection nullity.
- **OPEN:** whether a survivor, if any, is gauge, a pseudo-constraint or a
  physical propagating perturbation.

## Framing attack

This is not yet the Hessian *restricted as a quadratic form* to arbitrary
kinematic variations.  It asks the sharper first-order stationarity question:
which complete kinematic response vectors also lie in the action-selected
strong-equation graph?  A zero intersection closes this particular carrier
as a source of canonical evolution on the fixed slab.  A nonzero
intersection selects candidates but does not classify them physically.

The large raw scale coefficients near `lambda=1` can manufacture numerical
condition loss.  Rank decisions must use a frozen invertible column
equilibration, calibrated variant/radius errors, and an independent joined-
image formula.  Neither a raw condition number nor a desired continuum count
may select a threshold.

