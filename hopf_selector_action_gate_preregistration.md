# Preregistration: does the certified action contain the Hopf selector?

Date: 2026-08-10

## Complete question and hypotheses

The exact sixth-order invariant `S6(n)` selects the six Hopf axes only for a
potential with the sign `-g S6`, `g>0`.  This audit asks a narrower and
load-bearing question:

> Among the finite spectral triples and spectral moments already certified in
> the repository, is there a geometry-selected, `A5`-equivariant order
> parameter `n` and an internally generated sixth-order term with that sign?

The audit does not ask whether such a term can be added.  It asks whether it
is already produced by the theory's defined operators and admissible inner
fluctuations.

## Facts established before this audit

1. `verify_spectral_action.py` constructs the unfluctuated
   Kahler--Dirac operator and certifies only the finite moments
   `Tr(I)`, `Tr(D^2)` and `(1/2)Tr(D^4)`.  It defines neither a fluctuated
   `D_A`, a field `n`, nor a cutoff function with fixed coefficients.
2. On the free-cell arena, the canonical left `C[2I]` algebra has zero inner
   one-forms.  The right algebra has nonzero one-forms, but every enumerated
   real-structure candidate fails at least one required gate.  No physical
   fluctuated spectral action is licensed there.
3. The chamber witness `A=M2(C)+C^3` does pass the listed real finite-triple
   gates and has nonzero one-forms.  Its algebra type, cell support and
   colouring are explicitly labelled **STRUCTURAL**, not selected by the
   chamber geometry.
4. An exploratory calculation already found that the chamber block from the
   `(M2,C)` cell to the adjacent scalar cell has one-form complex dimension
   four and equals the image of `M2` multiplying the fixed Dirac block.  This
   is disclosed here; the result is not blind.

## Frozen exact tests

1. Independently rebuild the icosahedron, its 120 oriented chambers, the 60
   rotations of `A5`, the reflection, grading and chamber adjacency used by
   the B1 witness.
2. Reinsert the committed 60-entry B1 cell colouring.  Compute its exact
   stabilizer under all 60 chamber rotations.  Since its four cell capacities
   are distinct, no nontrivial permutation of cell labels is allowed when
   testing preservation of the represented central supports.
3. Rebuild the seven complex algebra-basis matrices.  On the noncommutative
   Dirac block, compute the exact span of

   `pi(a)[D,pi(b)]`, `a,b in A`.

   Compare it with `{(C tensor I_2)D_block : C in M2(C)}` on all matrix units.
4. If the coefficient module is `M2(C)`, verify its left `M2` action and the
   commuting right `M2` action.  As an `SU(2)` left module it is two copies of
   the fundamental doublet (its two columns), so a single doublet is not
   selected without an additional multiplicity vector.
5. Check the algebraic degree ceiling: when `D_A` depends linearly on any
   order parameter, `Tr(D_A^2)` and `Tr(D_A^4)` have polynomial degree at
   most two and four.  A nonzero homogeneous sixth-order anisotropy requires
   at least a `D_A^6` moment or a non-polynomial cutoff whose sixth coefficient
   has been specified.

## Decision boundary

The existing theory supplies the Hopf selector only if all of the following
are present simultaneously:

1. a valid real spectral triple with nonzero fluctuations;
2. a geometry-selected `A5`-equivariant three-vector or a uniquely selected
   doublet whose Hopf bilinear defines it;
3. a defined action containing a nonzero sixth-order invariant;
4. the sign that makes the six `C10` axes minima rather than the ten `C6`
   axes.

Failure of any item is a **DERIVED NEGATIVE for the current certified
construction**, not a theorem forbidding future dynamics.

In particular, a trivial `A5` stabilizer of the B1 colouring would show that
its doublets arise only after an embedding that has already broken the whole
icosahedral symmetry.  Such a witness cannot be used as evidence that the
geometry dynamically chooses one of six Hopf vacua.

The audit must keep separate:

- **DERIVED:** matrix ranks, module structure, stabilizer and degree bounds;
- **STRUCTURAL:** identifying a chamber `M2` doublet with a Higgs field;
- **OPEN:** a symmetry-preserving finite triple, the `D^6` coefficient and
  its sign, chirality locking, and continuum survival of the anisotropy.
