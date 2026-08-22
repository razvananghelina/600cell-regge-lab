# Prior-art gate: second finite-height full-boundary tangent and scale lift

Date: 2026-08-22

Status: **COMPLETED BEFORE CONSTRUCTING OR INSPECTING THE SECOND-SLAB
PRE-LEGENDRE MATRIX OR TANGENT MAP.**

## Question

For the fixed zero-`Lambda` Regge-plus-dust 600-cell history selected from
the representative incoming state `v=3/2`, construct the complete linearized
canonical map on the accepted second slab (branch B), express it in the same
physical canonical units as the first-slab map, and compose the two.

The immediate issue is not a desired spectrum.  It is whether scale
normalization at the intermediate slice can be undone without an arbitrary
canonical-variable choice.

## Primary literature checked

1. Dittrich and Hoehn, [Canonical simplicial
   gravity](https://arxiv.org/abs/1108.1974), derives discrete canonical
   evolution from Hamilton's principal function and emphasizes that adjacent
   regions compose by extremizing over their common boundary data.  This is
   the governing framework for the pre/post Legendre maps used here.
2. Dittrich and Hoehn, [From covariant to canonical formulations of discrete
   gravity](https://arxiv.org/abs/0912.1817), derives the canonical formalism
   directly from the discrete action and relates the linearized Hessian to
   constraints and gauge directions.
3. Bahr and Dittrich, [(Broken) Gauge Symmetries and Constraints in Regge
   Calculus](https://arxiv.org/abs/0905.1670), shows that curved Regge
   backgrounds generically replace exact constraints by pseudo-constraints.
   A regular finite-height Legendre map is therefore not by itself evidence
   that every boundary direction is a physical graviton.
4. Marsden and West, [Discrete mechanics and variational
   integrators](https://doi.org/10.1017/S096249290100006X), gives the standard
   discrete variational and generating-function origin of symplectic maps and
   their composition.
5. Grain and Vennin, [Canonical transformations and squeezing formalism in
   cosmology](https://arxiv.org/abs/1910.01916), explicitly documents that
   different time-dependent canonical variables can change finite-time mode
   and vacuum interpretations.  This is a direct warning against assigning
   physics to an isolated transfer-matrix spectrum before fixing the
   canonical variables geometrically.

The search also covered `Regge tangent map scale covariance`, `tent move
linearized perturbations`, `discrete Lagrangian cotangent lift scaling`, and
`time-dependent canonical cosmological perturbations`.  No source located in
this audit computes the present coefficient-level two-slab 600-cell dust map.
That search absence is not a novelty proof; external novelty remains **OPEN**.

## Repository facts already established

- **DERIVED EXACT / ADVERSARIALLY CORROBORATED:** under a global length
  scaling `L -> alpha L` with geometrized mass `M -> alpha M`, the complete
  action obeys `S -> alpha^2 S`.
- **DERIVED COMPUTATIONAL, MECHANICALLY DIFFERENTLY REPLICATED:** the first
  positive-height slab has a regular, symplectic 1440-dimensional canonical
  tangent, with no resolved even/odd staircase dependence.
- **DERIVED COMPUTATIONAL:** the accepted homogeneous branch B supplies a
  physical second slab and a unique physical third continuation.
- **OPEN:** physical constraint reduction, local propagating modes and any
  limiting speed.

## Framing correction

The proposed phrase "unique canonical co-moving trivialization selected by
scale covariance" is too strong.

Let `x_e=log|q_e|` and let `p_e=dS/dx_e`.  Under a global length scaling by
`alpha`,

```text
delta x -> delta x,
delta p -> alpha^2 delta p.
```

Thus

```text
D_alpha = diag(I, alpha^2 I)
```

satisfies `D_alpha^T Omega D_alpha=alpha^2 Omega`; it is
**conformal-symplectic**, not symplectic.  Calling it a canonical
transformation would be false.

Moreover scale covariance alone does not select a unique canonical frame.
Even after requiring H4 equivariance and identity at `alpha=1`, the families

```text
C_s(alpha)=diag(s(alpha) I, s(alpha)^(-1) I)
```

and the momentum shears

```text
H_b(alpha)=[[I,0],[b(alpha) I,I]]
```

are symplectic for arbitrary nonzero `s` and arbitrary `b`, with
`s(1)=1`, `b(1)=0`.  Therefore a single-step eigenvalue or singular value is
not made physical merely by calling one frame "co-moving".

## What scale covariance does select

It selects the physical-unit conjugacy of the complete tangent map.  If
`T_norm` is computed after dividing the entire second slab by its old length
scale `L_1`, then the direct physical-unit map must be

```text
T_phys = D_L1 T_norm D_L1^(-1),
D_L1   = diag(I, L_1^2 I).
```

Although `D_L1` itself is only conformal-symplectic, the conjugated map is
symplectic whenever `T_norm` is:

```text
T_phys^T Omega T_phys = Omega.
```

This is not a freely chosen trivialization.  It follows from the already
derived degree-two action scaling and can be falsified by assembling the
second-slab Hessian directly in physical units.  The shared intermediate
slice then uses literally the same physical edge labels and momenta, so the
two-slab differential is

```text
T_20 = T_2,phys T_1,phys.
```

Omitting `D_L1` can still produce a symplectic product, because it composes
two individually symplectic matrices in mismatched units.  Symplecticity is
therefore not a sufficient hostile control.  Direct physical-unit Hessian
agreement is mandatory.

## Acceptance and interpretation boundary

The next calculation advances only if:

1. the second normalized pre-Legendre system is regular;
2. its action-generated tangent is symplectic;
3. a separately assembled physical-unit Hessian equals `L_1^2` times the
   normalized Hessian within a frozen error ledger;
4. its direct physical-unit tangent equals
   `D_L1 T_norm D_L1^(-1)`;
5. both staircase schedules give the same physical map within the frozen
   classifier;
6. the four even/odd two-step compositions are mutually unresolved within
   propagated uncertainty.

Passing these gates proves a two-step linearized canonical response.  It does
not prove a physical mode spectrum.  Eigenvalues, continuum labels, a desired
wave equation, `c`, `G`, Planck units and particle data must remain unopened
until after this classification.

