# Adversarial nonhomogeneous direct-minor result

Date: 2026-08-20  
Status: **DERIVED COMPUTATIONAL; adversarially replicated; homogeneous sector OPEN**

## Frozen execution

The audit protocol was preregistered in `ad03ede`; the registered verifier source
was committed in `750f6ca`.

- verifier SHA-256:
  `459ffe5c7f6ac538c7682796b7bdffea891737ac2d727835fa8e00035af88bba`;
- artifact SHA-256:
  `ecf02fd76b0c1d4d95cd206c639a027400c2053bdb1850018d57ff2721861db3`;
- targeted execution: `7/7`;
- outcome: `NONHOMOGENEOUS_DIRECT_MINOR_REPLICATED`.

No full suite was run.

## Mechanically different certificate

The primary calculation used interval determinants of full Gram matrices. This
audit selected square row minors by deterministic pivoted QR and evaluated their
determinants directly with conservative Acb entry balls.

- The exact padded-identity positive control excludes zero.
- The exact duplicated-column negative control has determinant zero.
- The primary calculation reproduces `17/17` and rewrites a byte-identical frozen
  artifact.
- There are 12 nonhomogeneous parity-sector instances and 24 D/K matrices.
- All `48/48` cross-precision certificates exclude zero: P100-selected rows
  certify P160, and P160-selected rows certify P100.
- The frozen source/target reversal and pole-identity deletion falsifiers remain
  active.

Therefore the canonical scale--strut carrier/action intersection is **DERIVED
COMPUTATIONAL to be zero in every nonhomogeneous sector**, under the complete
frozen hypotheses of the primary construction.

## Boundary of the result

The audit intentionally does not retest the homogeneous block. Its no-refit
candidate misses the frozen exact-zero gate, so that block remains **OPEN**. No
classification as gauge, scale evolution, or physical propagation follows.

## Post-result primary-literature gate

The search was repeated using the technical terms learned in the calculation.

- Dittrich and Hoehn develop canonical simplicial gravity and action-generated
  Pachner evolution: [arXiv:1108.1974](https://arxiv.org/abs/1108.1974).
- Hoehn counts gauge and lattice-graviton degrees of freedom in linearized Regge
  calculus on flat backgrounds: [arXiv:1411.5672](https://arxiv.org/abs/1411.5672).
- Bahr and Dittrich show that curvature generically breaks exact discrete gauge
  symmetries into pseudo-constraints:
  [arXiv:0905.1670](https://arxiv.org/abs/0905.1670); Dittrich and Hoehn give the
  covariant-to-canonical construction in
  [arXiv:0912.1817](https://arxiv.org/abs/0912.1817).
- De Felice and Fabri evolve a generalized 600-cell by a Sorkin scheme:
  [arXiv:gr-qc/0106077](https://arxiv.org/abs/gr-qc/0106077).
- Liu and Williams study regular and perturbed closed Regge lattice universes:
  [arXiv:1502.03000](https://arxiv.org/abs/1502.03000).

**KNOWN:** nonhomogeneous simplicial evolution, lattice gravitons and perturbed
closed lattice universes exist in the literature. **OPEN external novelty:** this
search did not locate the same complete 600-cell scale--strut carrier, canonical
graph intersection, or its representation-sector direct-minor census. A web
search cannot prove absence, so no novelty claim is accepted.

## Consequence

This is a bounded negative about the present carrier/action compatibility, not a
derivation of gravity. The next exact question is confined to the single
homogeneous near-null direction. Until it is resolved, the route has no derived
tick, `c`, `G`, Planck scale, or mass formula.
