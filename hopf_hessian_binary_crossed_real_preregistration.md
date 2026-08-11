# Preregistration: ineffective binary crossed product and KO6 parity

Date: 2026-08-11

Parent effective crossed-product real no-go: `d84ea6d`.

## Why this arena must be checked

The five-point action factors through `A5`, but the theory's original group
is the binary icosahedral group `2I`.  Crossing by the effective quotient and
crossing by the full group with its ineffective central kernel are different
groupoid algebras.  The latter retains binary stabilizer data and must not be
silently identified with the former.

Fix

```text
B_bin=R(P) crossed_product 2I.
```

The stabilizer of a point is the preimage of `A4`, expected to be the binary
tetrahedral group `2T` of order 24.  This expected identification and every
real Wedderburn block must be verified from the exact binary group.

## STEP 1: exact stabilizer characters

1. reconstruct all 120 quaternion elements and the multiplication table;
2. identify the exact 24-element stabilizer/preimage `2T` and its normal
   `Q8` subgroup;
3. verify `2T/Q8=C3`;
4. construct the three one-dimensional characters;
5. construct the defining two-dimensional `SU(2)` character and its two
   `C3` twists;
6. construct the three-dimensional rotation character;
7. verify all seven complex characters are orthonormal and exhaust order 24;
8. compute every Frobenius--Schur indicator exactly.

The anticipated degrees

```text
1,1,1,2,2,2,3
```

are tests, not assumptions.

## STEP 2: real algebra and Poincare parity

Derive the complete real Wedderburn type of `R[2T]`, then use the transitive
groupoid identity

```text
B_bin ~= M5(R[2T]).
```

Check real dimension `5*120=600` and complexification block by block.  Count
the real simple summands and hence the free rank of `K0(B_bin)`.

Under the complete hypotheses

```text
KO6: J gamma=-gamma J,
metric-dimension-zero orientability,
strict Poincare duality over K0(B_bin),
```

test the generic antisymmetric intersection form at that rank.  Use a
symbolic rank/determinant argument, not a bounded multiplicity search.

## Decision boundary

- Even real `K0` rank: binary isotropy evades the previous parity no-go and
  the real Krajewski arena must be classified.
- Odd real `K0` rank: record a **DERIVED FULL-ARENA NO-GO** for the canonical
  ineffective binary crossed product, arbitrary multiplicities and `D`.

## Scope

This audit compares two canonical groupoids; it does not claim the
ineffective one is physically preferable.  A negative closes this natural
binary alternative under strict KO6 Poincare duality.  Dropping the central
kernel, Poincare duality or KO6 changes a stated hypothesis.

No Hessian or selector target is permitted.
