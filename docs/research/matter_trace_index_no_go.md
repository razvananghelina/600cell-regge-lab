# Matter trace-index test for the gauge prefactors

Date: 2026-07-22

## Exact target

Use `T(fundamental SU(n))=1/2`. For a finite chiral matter module

`M = direct_sum_i n_i (R3_i,R2_i)_{Y_i}`,

the trace indices are

`T1 = sum_i n_i dim(R3_i) dim(R2_i) Y_i^2`,

`T2 = sum_i n_i dim(R3_i) T(R2_i)`,

`T3 = sum_i n_i dim(R2_i) T(R3_i)`.

The required prefactors are proportional to `(8,5,2)` exactly when

`5 T1 = 8 T2`, `2 T2 = 5 T3`, with `T3!=0`.

For rational charges and integral multiplicities these are exact rational
Diophantine trace conditions.

## Standard-Model benchmark

For one left-handed Weyl generation

`Q=(3,2)_(1/6)`, `u^c=(bar3,1)_(-2/3)`,
`d^c=(bar3,1)_(1/3)`, `L=(1,2)_(-1/2)`, `e^c=(1,1)_1`,

ordinary hypercharge gives

`(T1,T2,T3)=(10/3,2,2)`, with ratio `5:3:3`, not `8:5:2`.
A sterile right-handed neutrino changes none of the indices. Multiplying the
abelian index by the conventional GUT factor `3/5` gives `(2,2,2)`, also not
the target.

- **DERIVED:** the desired coefficients are not the trace indices of one
  standard fermion generation under either of these standard conventions.
- **DERIVED (conditional statement):** if the ordinary one-generation content
  is retained and the target is scaled to its `T3=2`, additional color-singlet
  content would have to contribute exactly `(Delta T1,Delta T2,Delta T3)`
  `=(14/3,3,0)`. This is only a necessary trace condition; no such additional
  content is derived here.

## Inventory of existing finite spaces

1. **120 vertices and the `2I` regular-representation sectors.**
   **DERIVED:** these carry the regular `2I` actions and their irrep
   decomposition. **OPEN:** no action of all three derived gauge Lie factors
   on this space has been constructed.
2. **Affine-E8 McKay nodes.** **DERIVED:** these are the nine irreducible `2I`
   representations with tensoring by the defining `2I` doublet. They are not
   representations of the derived color algebra. A McKay-node label does not
   supply Standard-Model quantum numbers. The common gauge action is **OPEN**.
3. **The nine `(a,b)` fermion slots.** **DERIVED:** the slots have mass/grading
   assignments and particle-name identifications. They are labels, not a
   constructed vector space with three gauge actions. Their gauge-module
   decomposition is **OPEN**.
4. **The 12 Hopf fiber amplitudes.** **DERIVED:** before choosing brackets they
   carry the base `A5` permutation decomposition and Hodge data. With the
   compact bracket, `1+3+8` is the gauge algebra itself and therefore carries
   its adjoint action. **DERIVED conditional on that bracket choice:** its
   indices are `(0,2,3)`: `U(1)` acts trivially
   in the adjoint, while `T(adj SU(2))=2` and `T(adj SU(3))=3`. This is a gauge
   module, not a matter module, and it fails the target.
5. **Other scalar, vector, and Hodge spectral sectors.** **DERIVED:** their
   spectral and `A5`/`2I` decompositions are known where stated. **OPEN:** no
   common `u(1)+su(2)+su(3)` action is defined on them.

## Decision

- **DERIVED (negative):** no unconditionally derived matter candidate has a
  common action, so none supplies trace indices. Conditional on the compact
  color-bracket choice, the only testable common module is the gauge adjoint,
  with index ratio `(0,2,3)`, not `(8,5,2)`.
- **OPEN:** every matter-labelled candidate lacks the common gauge action
  required even to define `(T1,T2,T3)`.
- **PATTERN:** `(8/15,1/3,2/15)` remains a numerical prefactor pattern. It is
  not produced by a derived matter trace.

The missing object is now precise: an explicitly constructed finite module
`M`, built from already justified discrete data, together with matrices for
the three gauge-factor actions and a derived rational `U(1)` generator. Only
then can its three exact trace indices be compared with `(8,5,2)`.
