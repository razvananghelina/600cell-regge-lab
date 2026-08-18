# Preregistration: real and Galois forms of the six-fibration crossed product

Date: 2026-08-11

Prior complex algebra audit: `d15e7fa` protocol and its committed result in
`hopf_label_crossed_product_verdict.md`.

## Complete arena

Let `F` be the six derived Hopf fibrations with their exact transitive `A5`
action.  Audit both canonical coefficient forms of the transformation-group
algebra:

```text
B_R=R(F) crossed_product A5,
B_Q=Q(F) crossed_product A5.
```

The stabilizer is expected to be `D5` of order ten.  This and every character
and Wedderburn block must be reconstructed from the exact action.

No Hessian, `Box` point or selector target may be used.

## STEP 1: exact stabilizer representation theory

1. construct the stabilizer and all its conjugacy classes;
2. identify its rotation `C5` and reflection coset;
3. construct all four irreducible characters exactly in `Q(phi)`,
   `phi^2=phi+1`;
4. verify orthonormality and the degree-square sum;
5. compute all Frobenius--Schur indicators;
6. compute the exact golden-Galois action on the character set.

## STEP 2: coefficient forms

Derive independently:

- the real Wedderburn type of `R[D5]` and `B_R`;
- the rational Wedderburn type of `Q[D5]` and `B_Q`;
- real and rational dimensions;
- real/complex scalar extensions block by block;
- the free `K0` rank at each coefficient form.

The expected real four-block type

```text
M6(R)+M6(R)+M12(R)+M12(R)
```

is a test, not an assumption.

## STEP 3: KO6 and exact Galois compatibility

For the real split algebra, test only the antisymmetric parity gate; an even
rank means parity alone passes, not that a spectral triple exists.

For exact golden-Galois descent, let `P` swap the two conjugate `M12` blocks
and fix the two `M6` blocks.  Test both possible grading behaviors on a
generic four-node antisymmetric intersection form:

```text
PQP^T=Q,
PQP^T=-Q.
```

Record whether either permits nonzero Pfaffian.  Relate this to the rational
simple-summand count without conflating rational and real spectral triples.

## Canonical carriers

Record, without target comparison:

1. the minimum faithful left-module dimension and the commutant dimension of
   its standard graded double;
2. the regular bimodule and full enveloping bimodule standard doubles;
3. the already-derived natural six-label representation and whether it is
   faithful/order-zero capable.

An arbitrary proper Krajewski support does not count as canonical.

## Decision boundary

- If a canonical carrier passes and exact Galois descent permits Poincare
  duality, commit the blind census before comparing the Hessian.
- If the real algebra passes parity only after splitting the Galois pair, but
  no canonical carrier survives and exact descent is obstructed, record
  **STRUCTURAL REAL EXISTENCE / GALOIS SELECTION NEGATIVE** and stop.
- A full no-go may be claimed only with every coefficient-form hypothesis
  stated explicitly.

## Scope

Abandoning golden-Galois descent may leave abstract real Krajewski diagrams;
that is a changed arithmetic hypothesis.  Conversely, rational odd-rank
parity is not automatically a theorem about every real split triple.  Keep
these claims separate.
