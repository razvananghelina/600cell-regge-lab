# Complete real scale--strut carrier: symbolic gap resolved

Date: 2026-08-20

## Outcome

**DERIVED EXACT; ADVERSARIALLY CORROBORATED; KINEMATIC ONLY.**  On the real
domain

```text
lambda != 1,
tau != 0,
(lambda-1)^2-3 tau^2 != 0,
```

the compatible two-cell geometry uniquely fixes the endpoint response

```text
A = 6 - 2 tau^2/(lambda-1)^2,
B = 2 + 2 tau^2/(lambda-1)^2,
C = -1/(lambda-1),
D = lambda/(lambda-1).
```

Together with the three complete finite 600-cell controls, exact rank proof
and precision audit, this accepts the curved `1560 x 240` scale--strut map as
the complete kinematic carrier in the frozen coordinate choice.

It does **not** establish stationarity, an admissible symplectic phase space,
gauge directions, propagation, a physical tick, `c`, `G`, Planck scales,
gravitons, particle masses or Standard-Model physics.

## Frozen provenance

- gap-resolution protocol: `d6c62a7`;
- registered verifier: `260e89b`;
- frozen first run: `3d43113`, `10/11`, control failure;
- correction protocol: `e9d62b6`;
- assertion-only correction: `4fa0295`;
- resolved artifact: `438dca3`;
- resolved artifact SHA-256:
  `ea2c52f0cd227516734defc509330e528b140f71bfd0f50e87036f3fa9832179`;
- final targeted result: `11/11`,
  `FULL_SCALE_STRUT_GAP_REAL_RESOLVED`.

No full suite was run.

## What closed the gap

The pivot-free condition

```text
c_i r_j-c_j r_i = 0  for every i<j
```

reproduced the disclosed Gröbner ideal exactly.  Fixed-frame gluing still
gave the unit ideal and `D -> D+1` left 96 nonzero wedge residuals.

The complete actual denominator set was

```text
tau,
lambda-1,
(lambda-1)^2-3 tau^2.
```

The connection norm was exactly

```text
3*((lambda-1)^2+3 tau^2),
```

which is strictly positive over the stated real domain.  There were no
additional denominator or rank factors.  Direct exact rebuilds at
`(lambda,tau)=(-1,2)` and `(3,2)` showed that the former linear factors
`lambda+tau-1` and `lambda-tau-1` retain full local ranks, a rank-one
connection, and the disclosed ideal.  They were pivot artifacts, not real
exceptional strata.

The first corrected execution initially failed because the verifier wrongly
required `tau` to divide the connection norm.  That software failure and its
artifact are preserved separately.  The committed protocol never required
that false assertion; the correction changed only this control.

## Post-result primary-literature check

- Barrett et al., [*A Parallelizable Implicit Evolution Scheme for Regge
  Calculus*](https://arxiv.org/abs/gr-qc/9411008), already apply local vertex
  evolution to a 600-cell Friedmann cosmology.
- De Felice and Fabri, [*The Friedmann universe of dust by Regge Calculus:
  study of its ending point*](https://arxiv.org/abs/gr-qc/0009093), and their
  [generalized 600-cell evolution](https://arxiv.org/abs/gr-qc/0106077),
  establish that 600-cell Regge cosmology and less symmetry-restricted
  evolution are prior art.
- Dittrich and Hoehn, [*Canonical simplicial
  gravity*](https://arxiv.org/abs/1108.1974), make Hamilton's principal
  function/action—not kinematic compatibility—the generator of discrete
  canonical evolution.
- Hoehn, [*Canonical linearized Regge Calculus: counting lattice gravitons
  with Pachner moves*](https://arxiv.org/abs/1411.5672), identifies physical
  lattice modes through the action Hessian and emphasizes the special role
  of flat-background vertex-displacement symmetry.
- Tsuda and Fujiwara, [*Oscillating 4-Polytopal Universe in Regge
  Calculus*](https://arxiv.org/abs/2011.04120), derive symmetry-reduced
  4-polytopal FLRW dynamics and refinement limits.

The search found no primary source giving this exact endpoint-supported,
nonhomogeneous `120+120` response on a Lorentzian tetrahedral 600-cell slab.
Search absence is not a novelty proof.  **External novelty remains OPEN.**

## Framing attack and next gate

Calling this result “dynamics” would be false.  The map says which squared
length variations glue geometrically; it does not say which variations obey
the Regge equations or carry nonzero canonical momentum.

The next legitimate mission is therefore target-blind pullback of the frozen
Regge action derivative/Hessian to this accepted 240-dimensional carrier.
Before continuum labels are loaded it must record exact or controlled ranks,
nullities, symmetries and sector multiplicities, with corruption controls.
Only that calculation can distinguish gauge, constraint and genuinely
propagating directions.

