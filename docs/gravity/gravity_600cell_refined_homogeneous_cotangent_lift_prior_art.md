# Prior-art gate: homogeneous coarse-to-refined cotangent lift

Date: 2026-08-21

Status: completed before constructing the rank certificate.

## Exact object and complete hypotheses

Let the coarse spatial carrier be the regular 600-cell boundary and let

```text
K0 = P(sd K_600)
```

be the already certified projected barycentric refinement.  Restrict both
configuration spaces to their homothetic `H4`-invariant directions.  The
coarse boundary then has one spatial edge orbit, whereas `K0` has the six
rank-pair edge orbits

```text
01, 02, 03, 12, 13, 23.
```

Use logarithmic squared-edge coordinates.  If `s` is the common logarithmic
length scale, the refined configuration embedding has tangent

```text
h = dq_refined/ds = (2,2,2,2,2,2).
```

The question is whether the following data select a unique lift of a coarse
homogeneous momentum to the six refined orbit momenta:

1. the projected barycentric geometry;
2. full spatial `H4` invariance;
3. preservation of the canonical one-form, equivalently the cotangent
   pullback along the homothetic configuration embedding.

Both orbit-total and per-edge momentum conventions must be checked.  The two
are related by the already certified positive edge-population diagonal and
cannot change a rank or nullity.

This gate does not use an action Hessian, select an inner product on edge
space, solve a refined Regge slab, import a desired tick, or compare a mode
spectrum.

## KNOWN

- In discrete variational mechanics, pre- and post-momenta are obtained from
  the discrete action and the resulting evolution preserves the appropriate
  symplectic structure.  The cotangent construction supplies a pullback of
  covectors along a configuration map; it does not by itself supply a unique
  inverse when that map is an embedding.  See Marsden and West,
  *Discrete mechanics and variational integrators*, especially Sections
  1.4--1.5, DOI
  [10.1017/S096249290100006X](https://doi.org/10.1017/S096249290100006X).
- In canonical simplicial gravity, Hamilton's principal function is the
  action after internal variables have been eliminated; its boundary
  derivatives generate the canonical data.  See Dittrich and Hoehn,
  *Canonical simplicial gravity*, Sections 2 and 6,
  [arXiv:1108.1974](https://arxiv.org/abs/1108.1974).
- A refinement does not automatically define a dynamically faithful
  coarse/fine map.  Improved or perfect Regge actions obtain such relations
  by dynamical coarse graining.  See Bahr and Dittrich,
  *Improved and Perfect Actions in Discrete Gravity*,
  [arXiv:0907.4323](https://arxiv.org/abs/0907.4323), and Bahr, Dittrich and
  He, *Coarse graining free theories with gauge symmetries: the linearized
  case*, [arXiv:1011.3667](https://arxiv.org/abs/1011.3667).

These sources establish the framework, not the six-orbit rank result below.

## Repository controls

The refined feasibility census has already certified

```text
edge populations = (1440,3600,2400,3600,3600,2400),
six distinct H4 rank-pair edge orbits,
spatial edge count = 17040.
```

The exact refined rank geometry also certifies that all six squared edge
lengths are positive.  No schedule or Lorentzian slab is needed for this
boundary-only question.

## Framing attack

The canonical condition naturally runs from refined covectors to the coarse
homogeneous covector:

```text
p_s = h^T P = 2 sum_i P_i
```

for orbit-total momenta `P_i`.  Reversing this arrow is a lifting problem.
If the displayed functional has rank one, the refined momenta compatible
with a fixed `p_s` form an affine space of dimension five.  `H4` cannot
remove this freedom because each rank-pair orbit is separately invariant.

Choosing the Euclidean minimum-norm lift, population weighting, edge-length
weighting, a dual-volume weighting or an action supermetric can each select a
point, but those are additional structures.  They must be derived and named;
none follows from symplecticity alone.

Therefore the admissible negative conclusion is deliberately narrow:

> geometry, `H4` invariance and the canonical pairing alone do not select a
> unique homogeneous coarse-to-refined momentum lift.

It is forbidden to promote this to a theorem that no canonical transport can
exist after an action, a perfect-action coarse graining or an independently
derived supermetric is supplied.

## CONTROL / OPEN

- **CONTROL:** reproduce all six populations, their sum `17040`, and positive
  refined edge squares from the frozen exact rank formula.
- **CONTROL:** a one-orbit synthetic refinement has a unique momentum lift.
- **OPEN before calculation:** rank and nullity of the actual orbit-total and
  per-edge pullbacks.
- **OPEN:** whether the refined on-shell action supplies the missing five
  momentum components.
- **OPEN:** a refined physical phase space, tensor modes, dispersion, `c`,
  `G` and Planck units.

## Next admissible calculation

Preregister an exact symbolic verifier of the two pullback conventions and
the synthetic control.  If the nullity is five, direct reuse of the coarse
tick as refined initial data is structurally underdetermined.  The subsequent
physics route must obtain refined boundary momenta from a refined on-shell
action or an independently preregistered coarse-graining principle; it may
not choose one of the five directions by numerical convenience.

External novelty of the project-specific six-orbit statement is **OPEN**.
The literature search is not a novelty proof.
