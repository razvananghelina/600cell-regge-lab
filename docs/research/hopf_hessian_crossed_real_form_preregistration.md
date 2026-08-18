# Preregistration: real form and KO6 Poincare parity of the crossed product

Date: 2026-08-11

Parent complex-support census: `a21fd11`.

## Framing question

The complex crossed product has four summands

```text
M5+M5+M5+M15,
```

and exact coefficient conjugation exchanges the two summands induced from
the nontrivial `A4/V4=C3` characters.  The previous Pfaffian obstruction
assumed that a physical carrier respects this conjugation.

This audit asks whether that assumption is forced by the canonical real form
of the transformation-group algebra.  If so, the appropriate Poincare
pairing is indexed by real simple summands, not by four independently chosen
complex blocks.

## Complete algebraic arena

Fix the real transformation-group algebra

```text
B_R=R(P) crossed_product A5
   ~=M5(R[A4]),
```

where `P` is the derived five-point `A5` set.  Do not choose a different real
form after inspecting spectral-triple gates.

## STEP 1: exact real Wedderburn type

Using the exact constructed `A4` subgroup:

1. reconstruct its four conjugacy classes and quotient `A4/V4=C3`;
2. construct all four complex irreducible characters exactly in
   `Q(omega)`;
3. verify character orthonormality;
4. compute every Frobenius--Schur indicator

```text
nu(chi)=|A4|^-1 sum_h chi(h^2);
```

5. derive the real division-algebra type of every orbit of complex
   irreducibles;
6. derive `R[A4]`, `B_R`, their real dimensions and their complexifications.

Suggested values such as

```text
R[A4] ~= R+C+M3(R)
```

are expectations to test, not inputs to the verdict.

## STEP 2: KO6 intersection parity

Under the complete finite-triple hypotheses

```text
J gamma=-gamma J,
metric-dimension-zero orientability,
Poincare duality over K0(B_R),
```

verify the KO6 transpose sign of the intersection form and compute the
maximum possible rank for the number of real simple summands found in STEP 1.

No bounded multiplicity search is evidence for this claim.  Use the generic
antisymmetric matrix and its determinant/rank identity.

## Decision boundary

- If the real form has an even number of independent K0 generators, the
  parity argument does not close the arena; continue to real Krajewski types.
- If it has an odd number and the KO6 form is antisymmetric, record a
  **DERIVED FULL-ARENA POINCARE NO-GO** for the canonical real crossed product,
  arbitrary multiplicities and Dirac operators included.

## Scope

A no-go would depend on retaining the canonical real form, KO6 and Poincare
duality.  Treating the conjugate complex blocks as independent changes the
real algebra; dropping Poincare duality or changing KO dimension changes the
axioms.  State those exits plainly rather than calling them survivors.

No Hessian or selector target is permitted anywhere in this audit.
