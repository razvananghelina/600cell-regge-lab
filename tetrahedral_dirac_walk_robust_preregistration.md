# Preregistration: robust doubled tetrahedral walk on the H4 chamber graph

Date: 2026-08-11

## Provenance and motivation

The direct four-amplitude transplant was closed in commit `765c775`: natural
barycentric chamber orientation makes its causal stage two-to-one.

Appendix B of Nzongani *et al.*, *Dirac quantum walk on tetrahedra* (2024),
contains a distinct, previously published construction.  It doubles the
local carrier from four to eight amplitudes and factors the shift into three
swap stages.  The authors introduce it to make the walk robust under missing
links and irregular graphs.

This route is not invented or adjusted in response to our numerical output.
Its equations are frozen below before implementation.

## Complete operator class being tested

Use the already derived 14,400-chamber, four-coloured (H_4) graph.  Let
(n_2(k)) and (n_3(k)) be its intrinsic colour-2 and colour-3 involutions.
Give each chamber eight components (0,\ldots,7), for a total carrier
dimension

\[
8\times14{,}400=115{,}200.
\]

Implement Appendix B's three stages literally as output-to-input maps:

\[
S_0=(5,2,1,6,4,0,3,7),
\]

\[
S_1=(1,0,3,2,4,(n_2,6),(n_3,5),7),
\]

\[
S_2=(0,5,6,3,4,1,2,7).
\]

The macro shift is the synchronous global composition

\[
S=S_2S_1S_0.
\]

No free coefficient, coin angle, mass or phenomenological target is allowed.

## Exact tests

1. Reconstruct the same (H_4) chamber graph independently enough to check
   counts, colour involutions and the (3,3,5) Coxeter orders.
2. Check each of (S_0,S_1,S_2) is a permutation of all 115,200 states.
3. Check the macro composition is a permutation and crosses at most one dual
   chamber edge per macro tick.
4. Compute the macro map directly from the three global stages.
5. Separately transcribe the full-update column printed immediately after the
   three stage definitions in Appendix B and compare it with the literal
   composition.  Do not silently choose whichever version is unitary.
6. Record whether components 4--7 return to their initial chamber/component
   after a complete macro step and whether they are nevertheless used during
   the factorization.
7. Test independence from the global exchange of chamber orientation; the
   robust shift equations contain no handedness branch.

## Decision boundaries

- **DERIVED ROBUST SHIFT BRIDGE:** all three literal stages and their global
  composition are permutations on the (H_4) carrier.
- **DERIVED NEGATIVE:** any literal stage or composition is not bijective.
- **DOCUMENTED FORMULA DISCREPANCY:** the literal composition is unitary but
  differs from the paper's displayed expanded update.  In that case retain
  the operator definitions as the mathematically coherent construction and
  report the discrepancy explicitly; do not claim the printed expansion was
  verified.

## Scope boundary

Even a positive result establishes only a local unitary permutation
scaffold.  It does **not** show:

- equality with the Whitney Kähler--Dirac operator;
- a continuum Dirac limit on the curved 600-cell or its refinements;
- geometric selection of the Pauli-direction coin;
- a mass value, inertia, (c) in SI units or Planck scales.

The paper's flat-lattice continuum proof remains external and cannot be
ported by analogy.
