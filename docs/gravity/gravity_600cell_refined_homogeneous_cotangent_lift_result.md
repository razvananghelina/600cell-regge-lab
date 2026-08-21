# Consolidated result: coarse homogeneous momentum does not lift uniquely

Date: 2026-08-21

## Verdict

> **DERIVED EXACT / STRUCTURAL, adversarially corroborated:** on the
> homothetic sector of `K0=P(sd K_600)`, projected barycentric geometry, full
> spatial `H4` invariance and preservation of the canonical pairing determine
> only the pullback from six refined orbit momenta to one coarse homogeneous
> momentum.  The inverse fiber has exactly five free parameters.

This closes direct symplectic interpolation of the already accepted coarse
tick as a unique source of refined initial momentum.  It does not close an
action-selected or perfect-action coarse/fine transport.

## Provenance and reproduction

| stage | commit |
|---|---|
| prior-art and framing gate | `3e188e7` |
| primary target-free protocol | `3aefd29` |
| primary verifier registered before execution | `26e6ccc` |
| first primary result | `f83cdfa` |
| adversarial protocol | `282f570` |
| adversarial verifier registered before execution | `cbf4b7c` |

Primary verifier and artifact:

```text
reproducible/verify_gravity_600cell_refined_homogeneous_cotangent_lift.py
reproducible/gravity_600cell_refined_homogeneous_cotangent_lift.json
SHA-256 93dd857bff3b406e86d41a8a4b05d6441cb0e3e1c11e4f53d098555b1218924b
```

Adversarial verifier and artifact:

```text
reproducible/verify_gravity_600cell_refined_homogeneous_cotangent_lift_adversarial.py
reproducible/gravity_600cell_refined_homogeneous_cotangent_lift_adversarial.json
SHA-256 5489cd677c8414f69393209a95b24b5738d8fad9bbfb16462d7e4e25ae9f8a23
```

Both primary targeted runs were byte-identical and passed `12/12`.  Both
adversarial runs were byte-identical and passed `14/14`.  The static registry
audit reports `370` registered entries, `370` distinct names, zero duplicates,
zero unregistered verifiers, zero stale registrations and two reasoned
deliberate exclusions.  No full suite or nested `H4` census was run.

## Exact primary calculation

In logarithmic squared-edge coordinates the refined homothetic tangent is

```text
h=(2,2,2,2,2,2).
```

For orbit-total momenta,

```text
p_s = 2(P_01+P_02+P_03+P_12+P_13+P_23).
```

The exact rational pullback has

```text
rank = 1, nullity = 5.
```

For common per-edge momenta, the pullback row is

```text
2(1440,3600,2400,3600,3600,2400),
```

again with rank one and nullity five.  The population diagonal is invertible,
with determinant

```text
386983526400000000000,
```

so the two conventions are exactly equivalent.  Coordinate rescaling from
log squared length to log length and reversal of orbit order preserve the
result.  A synthetic one-orbit refinement has nullity zero and supplies the
required unique-lift positive control.

## Mechanically independent audit

The adversarial verifier imports no primary function and calls no rank,
nullspace, SVD or pseudoinverse routine.  For unit coarse momentum it builds
the six explicit lifts

```text
L_i=e_i/2.
```

Relative to `L_5`, the first five coordinates of the five differences form a
diagonal minor with exact determinant

```text
1/32.
```

The independently constructed per-edge differences have nonzero determinant

```text
1/5159780352000000000.
```

Z3 rational arithmetic finds two distinct six-orbit lifts with the same
coarse momentum (`SAT`) and proves that two distinct one-orbit lifts are
impossible (`UNSAT`).  The freedom persists at coarse momentum `-7/3` and
under orbit reversal; a corrupted `e_0/3` control correctly fails.

An intentionally stronger full permutation condition forces

```text
(1/12,1/12,1/12,1/12,1/12,1/12)
```

uniquely.  That is not an allowed repair: the actual barycentric colour
classes have distinct sizes `(120,720,1200,600)`, and spatial `H4` preserves
face dimension.  It fixes the six rank-pair labels rather than permuting them
as an `S6`.

## What is closed and what remains open

- **DERIVED EXACT / STRUCTURAL:** the declared geometry and symplectic pairing
  leave a five-dimensional affine freedom in refined homogeneous momentum.
- **CLOSED:** claiming that the coarse accepted tick uniquely initializes the
  refined `H4` slab by canonical interpolation alone.
- **FORBIDDEN WITHOUT A NEW DERIVATION:** Euclidean minimum norm, population,
  edge-length or dual-volume weighting chosen merely to pick one lift.
- **OPEN:** a lift selected by a refined on-shell action, an independently
  derived kinetic supermetric or a perfect-action coarse graining.
- **OPEN:** a refined physical phase space, tensor modes, dispersion, `c`,
  `G`, Planck scales and particle masses.

The result also sharpens the unequal-boundary framing.  Conserved dust mass
plus internal stationarity does not choose the outgoing boundary geometry,
and the coarse incoming momentum cannot be imported uniquely.  A refined
dynamic solve must include an action-derived refined initial-value condition;
otherwise its five hidden momentum choices are fitting.

## Post-result primary-source audit

Marsden and West's variational-integrator framework supplies the canonical
one-form and cotangent pullback, while Dittrich and Hoehn derive simplicial
pre/post momenta from Hamilton's principal function.  These are the relevant
known frameworks, not proofs of a reverse lift:

- Marsden and West, *Discrete mechanics and variational integrators*,
  DOI [10.1017/S096249290100006X](https://doi.org/10.1017/S096249290100006X);
- Dittrich and Hoehn, *Canonical simplicial gravity*,
  [arXiv:1108.1974](https://arxiv.org/abs/1108.1974).

The refined search again points to dynamical coarse graining when a map
between resolutions is required:

- Bahr and Dittrich, *Improved and Perfect Actions in Discrete Gravity*,
  [arXiv:0907.4323](https://arxiv.org/abs/0907.4323);
- Bahr, Dittrich and He, *Coarse graining free theories with gauge
  symmetries: the linearized case*,
  [arXiv:1011.3667](https://arxiv.org/abs/1011.3667).

No located primary source computes this exact six-orbit 600-cell pullback.
Search failure is not proof of novelty; external novelty remains **OPEN**.

## Next admissible physics gate

Do not solve an unequal-boundary refined slab against an arbitrarily lifted
coarse momentum.  The next construction must derive refined canonical data
from the refined action itself—for example through a preregistered
time-reflection-symmetric two-slab boundary-value problem—or derive a genuine
perfect-action coarse/fine map.  Its complete boundary conditions and matter
variables must be fixed before a root is inspected.

The preserved nested static census remains deferred at checkpoint `2/12`.
Its earlier complete no-root run is user-witnessed but not a certified
artifact; resume it only if a later conclusion depends on exclusion of that
equal-boundary branch.
