# Prior-art gate: precision resolution of the complete-carrier intersection

Date: 2026-08-20

## Disclosed target and exact scope

The binary64 census is frozen and numerically open.  Its post-result pattern
is now disclosed:

```text
six non-homogeneous minimal sectors:
    weaker scaled singular values around 1e-7, operational nullity zero;

homogeneous trivial sector:
    one reduced singular value around 3.9e-16,
    four further values around 1e-7,
    reduced/joined binary nullities disagree.
```

This mission is therefore target-disclosed.  It must try to refute both
tentative readings:

1. every non-homogeneous sector has zero intersection;
2. the homogeneous sector has exactly one intersection direction.

All hypotheses of the first census remain fixed: accepted curved
Regge--dust slab, both staircase parities, real nondegenerate scale--strut
formula, exact rank-240 carrier, frozen action Hessian, the same 65 response
orbits and seven deterministic minimal `2T` sectors.  No continuum
scalar/vector/tensor label is loaded.

## Primary-literature boundary

- Dittrich, Freidel and Speziale,
  [arXiv:0707.4513](https://arxiv.org/abs/0707.4513), relate exact Regge
  Hessian zero modes to remnant diffeomorphism symmetry and distinguish them
  from the continuum graviton propagator.
- Bahr and Dittrich,
  [arXiv:0905.1670](https://arxiv.org/abs/0905.1670), show that curved
  backgrounds generically turn exact discrete gauge constraints into
  pseudo-constraints.
- Hoehn, [arXiv:1411.5672](https://arxiv.org/abs/1411.5672), classifies
  lattice gravitons only after canonical/gauge analysis on flat
  backgrounds.
- Dittrich, Kaminski and Steinhaus,
  [arXiv:1404.5288](https://arxiv.org/abs/1404.5288), analyze geometric
  singularities and determinants of Regge Hessians.

Thus Hessian null modes, their gauge ambiguity and determinant tests are
**KNOWN**.  No located primary source resolves this exact 600-cell
scale--strut/action intersection.  Search absence does not establish
novelty; external novelty remains **OPEN**.

## Framing attack

The first error budget used a resolved difference between two *singular
spectra* of the global carrier as a conservative matrix perturbation proxy.
That was legitimate for forcing an OPEN result but cannot be reused as an
entrywise/matrix error theorem.  The resolver must construct the projected
carrier directly in multiprecision so that this proxy disappears rather
than merely lowering its coefficient.

A tiny midpoint singular value is not an exact kernel.  Full rank may be
certified by a nonzero interval Gram determinant.  Rank deficiency cannot be
certified merely because such a determinant contains zero.  For a proposed
one-dimensional kernel the calculation must separately establish rank at
least `n-1`, validate a lower-precision candidate at a higher precision
without refitting, and compare residual decay with the next singular scale
and all interval radii.  If those tests do not separate, the homogeneous
sector remains **OPEN**.

Even a resolved homogeneous direction would only be a first-order
action-compatible perturbation on one fixed slab.  It would not yet be a
gauge theorem, multi-tick evolution, emergent time, `c`, `G`, Planck units or
particle physics.

