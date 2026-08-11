# The six-fibration real algebra passes parity only after golden-Galois splitting

Date: 2026-08-11

Protocol commit: `beca527`.

Registered verifier:
`reproducible/verify_hopf_six_crossed_real_galois.py`.
Targeted exact result: `19/19`.

No Hessian or selector target was used.

## Headline

The six-fibration crossed product is the first canonical algebra in this
sequence whose split real form has even `K0` rank:

```text
B_R=M6(R)+M6(R)+M12(R)+M12(R).
```

Therefore KO6 antisymmetric parity alone does not kill it.

However, the two `M12` blocks are exchanged by the exact golden-Galois
automorphism.  Preserving the rational/Galois descent forces every KO6
intersection Pfaffian to vanish.  Moreover, none of the carriers furnished
functorially by the algebra survives the cheap real-triple gates.

The honest result is therefore:

- **STRUCTURAL REAL EXISTENCE:** even-rank split real algebra;
- **DERIVED GALOIS-DESCENT NO-GO:** no nondegenerate Galois-compatible KO6
  intersection form;
- **DERIVED CANONICAL-CARRIER NEGATIVE:** zero functorial carriers reach a
  Hessian comparison.

## Exact stabilizer theory

The exact `A5` action on the six fibrations is transitive with point
stabilizer `D5` of order ten.  Its conjugacy-class size/order data are

```text
sizes  =1,5,2,2,
orders =1,2,5,5.
```

The four exact characters in `Q(phi)`, `phi^2=phi+1`, have degrees

```text
1,1,2,2
```

and are orthonormal.  All Frobenius--Schur indicators equal `+1`, so every
complex irrep is of real type.  Golden Galois fixes the trivial and reflection
sign characters and exchanges the two doublets.

## Three coefficient forms

Over the reals,

```text
R[D5] = R+R+M2(R)+M2(R),
B_R   = M6(R)+M6(R)+M12(R)+M12(R).
```

The real dimension is `36+36+144+144=360`, and the real free `K0` rank is
four.

Over the rationals, the Galois pair remains one field block:

```text
Q[D5] = Q+Q+M2(Q(sqrt(5))),
B_Q   = M6(Q)+M6(Q)+M12(Q(sqrt(5))).
```

Its rational dimension is `36+36+288=360`, but it has three rational simple
summands.  Complexification gives the already-derived type

```text
M6(C)+M6(C)+M12(C)+M12(C).
```

These statements are compatible: the field block has two real/complex
embeddings, which become the two `M12` summands only after scalar splitting.

## KO6 and golden Galois

For four independent real nodes, the generic KO6 intersection form has

```text
Pf(Q)=Q01 Q23-Q02 Q13+Q03 Q12,
det(Q)=Pf(Q)^2,
```

so nondegeneracy is algebraically possible.

Exact Galois exchanges only the two `M12` nodes.  Its node permutation `P`
is a transposition with determinant `-1`, hence

```text
Pf(PQP^T)=-Pf(Q).
```

If Galois preserves grading, `PQP^T=Q` forces `Pf(Q)=0`.  If it reverses
grading, `PQP^T=-Q`, while `Pf(-Q)=Pf(Q)` in four dimensions, with the same
result.  The verifier solves both generic constraint spaces exactly.

Equivalently, before splitting the golden field the rational form has three
simple blocks, so its antisymmetric intersection pairing is necessarily
degenerate.

This is a no-go for exact Galois descent, not for every real split algebra
obtained after choosing one embedding of the golden field independently.

## Canonical carriers

The minimum faithful left module has dimension

```text
6+6+12+12=36.
```

Its standard graded double has dimension 72 but commutant dimension only 16.
A faithful 360-dimensional opposite algebra cannot fit there, so order zero
fails.

The regular and full-enveloping standard doubles represent every metric-zero
cycle with sheet profile `(1,1)`, while `gamma` has `(1,-1)`; both fail
orientability.

The natural six-label representation has image full `M6` and kernel
dimension 324.  Its same-branch odd double fails order zero.  Its `A5`
restriction is the constant line plus the physical five-dimensional module;
the full `M6` algebra does not preserve that physical summand separately.

Thus the ambient six-label construction does not provide a licensed
12-state home for the affine Hessian merely by adjoining its constant zero
mode.

## Why no Hessian comparison follows

The preregistered comparison gate required both an independently selected
carrier and exact arithmetic compatibility.  No canonical carrier survives,
and the only available exact Galois descent forbids Poincare duality.

One could abandon Galois descent and enumerate proper four-node real
Krajewski supports.  That would be an abstract existence search with a new
arithmetic choice, not a consequence of the current geometry.  Selecting a
support by its later Hessian behavior would be fitting.

## Status ledger

- **DERIVED:** the exact stabilizer is `D5`, with character degrees
  `1,1,2,2` and FS indicators all `+1`.
- **DERIVED:** `B_R=M6(R)^2+M12(R)^2`, real dimension 360 and even `K0` rank
  four.
- **STRUCTURAL OPENING:** KO6 parity alone permits a nondegenerate real split
  intersection form.
- **DERIVED:** `B_Q=M6(Q)^2+M12(Q(sqrt5))`, with three rational simple blocks.
- **DERIVED GALOIS NO-GO:** zero nondegenerate intersection forms preserve or
  reverse grading under exact golden conjugation.
- **DERIVED CANONICAL-CARRIER NEGATIVE:** the minimum faithful left double,
  regular double, enveloping double and natural six-label double all fail a
  necessary gate.
- **NO TARGET COMPARISON:** the Hessian was not inspected.
- **OPEN:** a geometry-derived mechanism that breaks golden Galois and
  selects a proper real bimodule before target comparison.

## Programme consequence

The six-fibration algebra is closer than the five-point algebras: its split
real type has the required even number of summands.  But the improvement is
not free.  It occurs exactly by splitting a conjugate golden pair, and the
current geometry supplies no reason to choose that split while discarding
Galois descent.  Until such a mechanism is derived, this is a structural
opening rather than a physical finite triple.
