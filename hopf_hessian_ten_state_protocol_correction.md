# Protocol correction: the five-dimensional A5 irrep admits a normalized C5 algebra

Date: 2026-08-11

Parent preregistration: `a0d7dd8`.

## The error caught before gate computation

The preregistration proposed the sheet-image list

```text
C I_5, M5(C)
```

and suggested that a five-support algebra would force the reducible pure
permutation module `1+4`.  That implication is false for a normalized
diagonal algebra.  Its unitary normalizer consists of **monomial** matrices,
not only permutation matrices; stabilizer phases can remove the fixed vector.

For `A5`, an index-five stabilizer is `A4`, whose abelianization is `C3`.
For either nontrivial character `chi:A4->C^*`,

```text
Ind_A4^A5(chi)
```

has dimension five and is the irreducible `W_5`.  Thus `W_5` normalizes a
maximal diagonal algebra `C^5` in an induced monomial basis.

This is a direct counterexample to the incomplete expected list.  The
ten-state no-go cannot be inferred from `dim End_A5(W_5)=1`.

## Corrected STEP 1

The sheet-image classification to verify is now

```text
C I_5,
C^5 in an A5-monomial basis,
M5(C).
```

The central-support argument still excludes centre sizes two, three and four:
`A5` has no nontrivial permutation action of those degrees, so the supports
would be individually invariant and contradict irreducibility.  Centre size
five is retained and must be constructed explicitly.

## Required exact construction

Before applying any spectral-triple gate:

1. construct an `A4` subgroup of the exact 60-element action;
2. construct both nontrivial `A4` character projectors;
3. verify that each projector has rank one;
4. take its `A5` orbit and verify exactly five mutually orthogonal rank-one
   projectors summing to `I_5`;
5. verify that the group normalizes their `C^5` span.

Only then may `C^5` embeddings be included in the joint-type and first-order
enumeration.

## Corrected first-order task

For a commutative five-point algebra, first order constrains a Dirac block to
a Krajewski rook pattern determined by the left and opposite label maps.  It
does not automatically force the whole block to zero.

Enumerate all relative permutations allowed by the two sheet embeddings and
the KO6 sheet-exchanging `J`.  For every case, compute the exact allowed
matrix-entry support and compare the **entire affine family**

```text
bI+cHhat_X, c!=0,
```

with that support.  A test in one preferred monomial basis is insufficient;
the conjugate character system and all relative label permutations must be
covered.

## Revised decision boundary

- If a `C^5` embedding supports the complete Hessian family, continue through
  nonzero forms, connectedness, orientability and Poincare duality.
- If every exact rook support misses some coefficient of the five-dimensional
  Hessian family, close only the fixed ten-state carrier.
- Full `M5` types may still fail order zero and scalar `C` may still have zero
  forms, but those facts cannot substitute for the missing `C^5` audit.

## Provenance

This correction is committed before construction of the monomial projectors
and before any first-order support result.  The original preregistration is
left in history so the mistake and its repair remain externally visible.
