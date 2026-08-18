# Preregistration: the canonical five-point crossed product

Date: 2026-08-11

Parent commutative no-go: `499f9ff`.

## Complete arena

Let `P` be the five primitive projectors of either exact monomial `C^5`
system on `W_5`.  The two systems have the same derived transitive
permutation action of `A5` on `P`.  Fix the transformation-group algebra

```text
B=C(P) crossed_product A5.
```

This is the next noncommutative algebra selected by the existing geometry.
No proper subalgebra, central summand or Krajewski support is to be selected
after inspecting the Hessian.

## STEP 1: target-blind algebra and carrier census

Without constructing or comparing `Hhat_X`:

1. reconstruct the exact five-point action and its stabilizer `H`;
2. compute the conjugacy classes and complex irreducible degrees of `H`;
3. derive the complete Wedderburn type of `B` from
   `B ~= M_5(C[H])`;
4. construct the covariant representations induced from every
   one-dimensional stabilizer character and compute their exact images and
   kernels;
5. identify which of those representations are the trivial five-point
   permutation module and the two monomial `W_5` systems;
6. compute the minimum dimension of a faithful left `B`-module;
7. audit the two carriers furnished without an extra Krajewski choice:
   the regular bimodule `L^2(B)` and the full enveloping bimodule
   `B tensor B^op`.

Record all results and commit them before any Hessian comparison.

## Canonical-carrier criterion

A carrier counts as canonical in STEP 1 only if it is functorially supplied
by `B` itself or by one of the already-derived stabilizer characters.  An
arbitrary subset of Wedderburn bimodule cells does not count merely because
it satisfies order zero.

For the regular bimodule, record its central cell support.  For the full
enveloping bimodule, record whether diagonal central cells occur.  Test the
KO6 metric-dimension-zero orientability sign equation on those supports.

## STEP 2: only after the blind commit

If a canonical carrier survives faithfulness, order zero and orientability,
then compare its first-order map space with the complete affine Hessian
family.  Record every induced-character branch and the hit fraction.

If no canonical carrier survives, stop without manufacturing a smaller
Krajewski graph.  The honest verdict is then:

```text
the crossed-product algebra is selected, but it does not select a viable
physical carrier.
```

## Expected values are tests, not assumptions

Groupoid theory suggests, but the verifier must establish, that `H=A4` and

```text
B ~= M5(C[A4]) ~= M5+M5+M5+M15.
```

It also suggests that the three five-dimensional covariant representations
select the three `M5` summands and that each monomial `W_5` image is full
`M5`.  A disagreement changes the route and must be reported rather than
patched.

## Framing boundaries

- A full `M5` image on `W_5` is not itself a physical success; order zero and
  faithfulness still matter.
- The existence of some abstract faithful bimodule is not canonicity.
- The regular and enveloping carriers may fail orientability.  That closes
  those carriers, not every abstract Krajewski diagram over `B`.
- A later choice among multiple central-cell supports must carry an explicit
  look-elsewhere count and an independent geometric selection rule.
