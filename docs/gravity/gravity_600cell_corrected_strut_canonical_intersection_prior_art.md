# Prior-art gate: corrected-strut/canonical-graph intersection

Date: 2026-08-19

## Exact object and hypotheses

Fix the accepted non-static Lorentzian 600-cell Regge--dust slab, separately
for both staircase parities.  In each of the seven minimal `2T` sectors of
dimension `d`, fix:

- the target-blind geometry-selected corrected strut graph
  `G : C^(5d) -> C^(65d)` from commit `dab941b`;
- the canonical strong-equation graph
  `C = (-A^-1 B,I)^T : C^(5d) -> C^(65d)` reconstructed from the frozen
  action Hessian.

The five pole positions and their coefficient order are fixed by edge
geometry.  Both graphs have the literal identity on those same pole rows.
Consequently, a vector belongs to both graph images iff it has one common
coefficient vector `x` and

```text
(G-C)x = 0.
```

The new object is the nullity and kernel of `G-C` in all
`2 parities x 7 sectors = 14` cases.  No physical mode label or expected
nullity is assumed.

## Primary-literature boundary

- Dittrich and Hoehn,
  [*From covariant to canonical formulations of discrete gravity*,
  arXiv:0912.1817](https://arxiv.org/abs/0912.1817), derive canonical
  constraints from the Regge action and explain their deformation into
  background-dependent pseudo-constraints on curved backgrounds.
- Dittrich and Hoehn,
  [*Canonical simplicial gravity*,
  arXiv:1108.1974](https://arxiv.org/abs/1108.1974), formulate discrete
  evolution through pre/post constraints and Hamilton's principal function.
- Hoehn,
  [*Canonical linearized Regge Calculus: counting lattice gravitons with
  Pachner moves*, arXiv:1411.5672](https://arxiv.org/abs/1411.5672), counts
  vertex-displacement and curvature-carrying modes on flat backgrounds.
- De Felice and Fabri,
  [arXiv:gr-qc/0009093](https://arxiv.org/abs/gr-qc/0009093) and
  [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077), study discrete
  dust evolution of the 600-cell but not this anisotropic graph intersection.
- Jercher and Steinhaus,
  [arXiv:2312.11639](https://arxiv.org/abs/2312.11639), study dynamical
  Lorentzian frustum heights and matter-supported branches in a different
  symmetry-reduced carrier.

These sources make action-generated constraint graphs and their possible
failure to coincide with kinematic vertex-displacement data **KNOWN**.  None
of the located primary sources provides this corrected 120-column graph or
its sectorwise intersection with this slab's Hessian graph.  Search absence
does not prove novelty; external novelty is **OPEN**.

## KNOWN / CONTROL / OPEN

- **DERIVED INPUT:** the corrected carrier is rank 120, target-blind and
  geometry-selected.
- **DERIVED INPUT, adversarially corroborated:** its complete sector images
  are not equal to the canonical graph images.
- **CONTROL:** literal equality of the pole blocks makes `ker(G-C)` exactly
  the graph intersection coefficient space; without this equality a stacked
  image-intersection calculation would be required instead.
- **CONTROL:** `G=C` must give nullity `5d`; a frozen injected full-column-rank
  difference must give nullity zero.
- **OPEN:** every actual nullity and kernel.
- **OPEN:** whether a surviving intersection direction is gauge,
  pseudo-constraint or physical; nullity alone cannot decide that.

## Framing attack

The preceding projector distance is the sine of the largest principal angle.
It proves non-equality but cannot determine intersection dimension.  Reading
zero intersection from a distance near one would be a mathematical error.
This mission exists specifically to close that gap.

Even a nonzero kernel would not show that arbitrary strut data are lapse
freedom.  It would select only the common subspace.  A zero kernel would be a
clean negative for pure-strut canonical freedom, but would not decide the
full 240-dimensional scale-plus-strut carrier.

