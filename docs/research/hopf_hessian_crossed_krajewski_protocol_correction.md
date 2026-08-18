# Protocol correction: binary support does not fix weighted Poincare rank

Date: 2026-08-11

Parent preregistration: `8b21830`.

## Error caught before the STEP 1 result commit

The preregistration claimed that binary multiplicities were exhaustive for
minimum Hilbert dimension because duplicating a cell did not change
intersection rank.  That statement is false in general.

For four nodes,

```text
Pf(Q)=Q01 Q23-Q02 Q13+Q03 Q12.
```

Two nonzero terms can cancel when every occupied edge has multiplicity one
and cease to cancel at unequal positive multiplicities.  Thus a support with
singular binary `Q` can acquire rank four after weighting.

The binary `3^6=729` census remains a complete topology census at unit
multiplicity.  It is not by itself a proof of the unrestricted minimum.

## Corrected exhaustive minimum protocol

The first blind run supplies a unit-multiplicity survivor of total Hilbert
dimension 300.  Use 300 only as an upper bound, not as the answer.

Now enumerate every orientable signed multiplicity assignment with total
dimension at most 300:

```text
sum_(i<j) 2 n_i n_j |q_ij| <= 300,
n=(5,5,5,15),
q_ij in Z,
q_ji=-q_ij.
```

For the three `5-5` pairs this implies `|q_ij|<=6`; for the three `5-15`
pairs it implies `|q_ij|<=2`.  The finite box has

```text
13^3 * 5^3 = 274625
```

assignments before the dimension cut.  Enumerate it exhaustively and apply
the same faithfulness, Pfaffian-rank, first-order legality and central
connectedness gates to the nonzero support.

Record the exact minimum dimension and every weighted minimizer.  This
repairs the minimum claim but does not turn the infinitely many unrestricted
positive multiplicities into a finite look-elsewhere denominator.

## Multiplicity-independent conjugation test

Let `P` exchange the two conjugate `M5_+` and `M5_-` nodes.  It is a single
transposition, so `det(P)=-1`.  For every antisymmetric four-by-four matrix,

```text
Pf(P Q P^T)=det(P) Pf(Q)=-Pf(Q).
```

Test both possible grading behaviors:

1. grading preserving: `PQP^T=Q`;
2. grading reversing: `PQP^T=-Q`.

In four dimensions `Pf(-Q)=Pf(Q)`.  Either behavior therefore forces
`Pf(Q)=0`, hence degeneracy.  If verified symbolically, this closes every
conjugation-compatible multiplicity assignment, not merely the binary
census.

## Revised decision boundary

- The raw binary survivor count measures support-level nonuniqueness only.
- The corrected weighted search establishes the true minimum up to the
  certified upper bound.
- If the Pfaffian conjugation obstruction holds, no survivor is canonical
  under the exact character conjugation.  A non-invariant graph remains an
  abstract **STRUCTURAL EXISTENCE** witness, not a selected carrier.

No Hessian or selector comparison is permitted in this correction.
