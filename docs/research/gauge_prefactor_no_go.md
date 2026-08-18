# Gauge-prefactor status: normalization no-go

Date: 2026-07-22

## Decision

The coefficients `(8/15, 1/3, 2/15)` remain **PATTERN**. They are not derived
by the presently defined 600-cell operators.

## What is derived

- **DERIVED:** `ker(Box_1) = rho_0 + 2 rho_5`, with a 12-dimensional
  nontrivial part.
- **DERIVED:** after passage to the Hopf base, its `A_5` content is
  `1 + 3 + 3' + 5`, and the dimension grouping `1 + 3 + (3' + 5) = 1 + 3 + 8`
  matches the dimensions of `u(1) + su(2) + su(3)`.
- **DERIVED:** the available Hodge quadratic forms are scalar on the entire
  12-dimensional gauge kernel: `C = 5 I`, `B = (16/5) I`, and
  `Delta_1 = (41/5) I` there.

## No-go statement

The present discrete data do not select the three Yang--Mills trace
normalizations.

1. The `1`, `3`, and `8` sectors are inequivalent symmetry sectors. An
   invariant positive quadratic form can therefore be independently rescaled
   on each sector. Symmetry alone leaves three positive coefficients (two
   ratios after removing one overall scale).
2. The operators actually computed on `ker(Box_1)` do not remove this freedom:
   their restrictions are scalar, so their traces only count dimensions.
3. The color bracket problem has now been classified exactly. There is a
   unique compact-simple `A_5`-equivariant bracket class on `3'+5`, namely
   `su(3)`, and its Killing form fixes the relative normalization of the
   `3'` and `5` blocks. Equivariance and Jacobi also allow split and degenerate
   classes. The canonical edge metric further removes every noncompact and
   nonreductive noncentral branch, as described below.
4. A finite heat trace or spectral moment of these scalar restrictions cannot
   distinguish the three required normalizations without inserting additional
   projectors and their weights. Such inserted weights are precisely the
   missing information.

Consequently, obtaining `(8/15, 1/3, 2/15)` from a trace currently requires a
choice of sector weights equivalent to assuming the desired normalization.
That is not a derivation.

## Minimal missing structure

The bracket part has been reduced by `a5_equivariant_bracket_theorem.md`:

- **DERIVED:** the unique compact-simple bracket class on `3'+5` is `su(3)`;
- **DERIVED:** the `3` sector has the cross-product `su(2)` bracket up to scale;
- **DERIVED:** the `1` sector is abelian;
- **DERIVED:** the color Killing form fixes the internal `3':5` relative norm.
- **DERIVED:** the canonical edge metric is `20 Frob` on `3'` and `16 Frob`
  on `5`; ad-invariance selects only the abelian, `so(3)+R^5` central, and
  compact `su(3)` branches. Its compact representative has parameters
  `(a,0,a,-4a/5)`.

What remains for color is the single **STRUCTURAL** compatibility axiom that
its bracket makes the color sector a center-free metric Lie algebra for the
canonical edge inner product. Trivial center removes the abelian and
`so(3)+R^5` branches; ad-invariance removes the split, semidirect, and
nilpotent branches. What remains for the prefactors is a normalized
representation of `u(1)+su(2)+su(3)` on a common finite matter space,
including:

- an embedding of each sector into endomorphisms of that matter space, and
- one common trace convention (including the `U(1)` charge normalization).

Only then are the three coefficients computable as trace indices. Until such
data are constructed independently, the fractions are **PATTERN**, not
**STRUCTURAL** and not **DERIVED**.

## Finite matter-space audit

The follow-up computation in `matter_trace_index_no_go.md` makes the required
trace condition exact. With `T(fundamental SU(n))=1/2`, the target is

`5 T1=8 T2`, `2 T2=5 T3`, `T3!=0`,

equivalently `(T1:T2:T3)=(8:5:2)`.

- **DERIVED:** one ordinary-hypercharge Standard Model Weyl generation has
  `(T1,T2,T3)=(10/3,2,2)`, hence ratio `5:3:3`, not `8:5:2`.
- **DERIVED conditional on the compact bracket choice:** the only repository
  space currently carrying explicit actions of all three factors is the gauge
  algebra in its adjoint module. Its indices are `(0,2,3)`, not the target.
- **OPEN:** the 120-vertex regular sectors, McKay nodes, `(a,b)` fermion slots,
  and the other spectral sectors have no derived common gauge action. Their
  three indices are therefore undefined, not failed numerical candidates.
- **DERIVED (negative):** no unconditionally derived matter module even has
  all three actions; the conditionally available adjoint module fails.

Thus the obstruction is no longer an unspecified search over existing finite
spaces. A new construction is required: explicit matrices for a common matter
module and a rational `U(1)` generator derived from the discrete data.

## Open question

- **OPEN:** derive rather than postulate the center-free metric-Lie-algebra
  compatibility axiom for color.
- **OPEN:** construct the common matter representation/trace package and test
  whether its indices equal `(8/15, 1/3, 2/15)`.
