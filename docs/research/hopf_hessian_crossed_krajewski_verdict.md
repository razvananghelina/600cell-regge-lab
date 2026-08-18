# Abstract crossed-product Krajewski diagrams exist, but geometry selects none

Date: 2026-08-11

Initial protocol: `8b21830`.

Protocol correction: `ed4b989`.

Target-blind enumeration commit: `a21fd11`, explicitly made with no Hessian
comparison.

Verifier: `reproducible/verify_hopf_hessian_crossed_krajewski_blind.py`.
Targeted exact result: `12/12`.

## Verdict

Two statements must not be conflated.

1. **STRUCTURAL EXISTENCE:** the algebra
   `M5+M5+M5+M15` admits many abstract central Krajewski supports satisfying
   the audited necessary combinatorial gates.
2. **DERIVED SELECTION NEGATIVE:** none is compatible with the exact
   conjugation exchanging the two nontrivial stabilizer-character blocks,
   and without that symmetry there are hundreds of attempts rather than a
   selected carrier.

Therefore no STEP 2 Hessian comparison was performed.  Choosing one of the
survivors because it later fits the selector would violate the preregistered
acceptance boundary.

## Look-elsewhere result

The complete unit-multiplicity support space has 729 elements.  The exact
necessary-gate ledger is

```text
Poincare rank four                 484
faithful                           636
nonzero first-order support        642
central connectedness              316
all gates                          256
```

The raw survivor fraction is

```text
256/729 = 35.1 percent.
```

This is not a rare hit.  Even after quotienting by the full abstract `S3`
permutation of equal `M5` summands and by grading reversal, 22 classes remain.
That quotient grants more symmetry than the geometry, because the trivial
stabilizer character is distinguished from the two irreducible `W5`
branches.

The true weighted minimum was audited after correcting an error in the
initial binary argument.  Exhausting every integer multiplicity assignment
of dimension at most the binary upper bound gives

```text
assignments in signed search box      274625
assignments after dimension cut           773
weighted survivors                         24
minimum total Hilbert dimension            300
minimum designs                             24.
```

Every minimizer is binary.  Thus the value 300 is proved, but it is attained
24 ways and supplies no uniqueness.

## Exact symmetry obstruction

Let `P` exchange the central nodes belonging to `chi` and `chibar`.  Since
`P` is a transposition, `det(P)=-1`.  For every alternating four-node
intersection form,

```text
Pf(PQP^T)=-Pf(Q).
```

If conjugation preserves grading, it requires `PQP^T=Q` and hence
`Pf(Q)=0`.  If it reverses grading, it requires `PQP^T=-Q`; but
`Pf(-Q)=Pf(Q)` in four dimensions, again giving `Pf(Q)=0`.

This proof covers arbitrary integer multiplicities and is independent of the
729-support search.  The exact symbolic verifier confirms both generic
constraint spaces.

Hence:

> **DERIVED CONJUGATION NO-GO.** Under either grading behavior, no
> nondegenerate Krajewski intersection form over the four crossed-product
> summands is compatible with the exact `chi<->chibar` conjugation.

The reason this conjugation is canonical has now been derived independently
in `hopf_hessian_crossed_real_form_nogo.md`.  Exact Frobenius--Schur indicators
give the real algebra

```text
M5(R)+M5(C)+M15(R),
```

whose complex-type middle block splits into `chi,chibar`.  Its real `K0` rank
is three, so strict KO6 Poincare duality fails directly by odd antisymmetric
rank.  The Pfaffian obstruction is the complexified shadow of this real-form
no-go, not an arbitrary extra symmetry demand.

## What the survivors do and do not show

The 256 binary survivors show only that dropping conjugation leaves abstract
central designs.  They have not passed:

- construction of an antiunitary on the complete matrix factors;
- a numerical self-adjoint Dirac operator;
- full first order for that operator;
- nonzero represented one-forms;
- connectedness beyond the centre;
- a derived coupling to the five-dimensional Hessian family.

The central census could falsify existence, but it could not certify these
later gates.  It is therefore not cited as evidence for a complete spectral
triple.

## Programme status

- **DERIVED ALGEBRA:** the five-point crossed product fixes
  `M5+M5+M5+M15`.
- **STRUCTURAL EXISTENCE:** `256/729` binary central supports pass necessary
  gates.
- **PATTERN/LOOK-ELSEWHERE NEGATIVE:** 24 minimum designs and 22 generously
  quotiented classes remain.
- **DERIVED CONJUGATION NO-GO:** zero nondegenerate designs preserve or
  reverse grading under exact character conjugation, arbitrary multiplicity.
- **DERIVED REAL-FORM NO-GO:** the canonical real crossed product has three
  simple summands and cannot carry a nondegenerate KO6 intersection form.
- **DERIVED TARGET-PROTOCOL STOP:** no independently selected carrier
  survives, so comparing the Hessian would be fitting.
- **STRUCTURAL SELECTOR RETAINED:** the exact affine fourth-moment Hopf
  selection remains true as operator algebra.
- **OPEN PHYSICS:** a different geometry-derived real algebra/carrier or a
  justified change of KO/orientability hypotheses.

## Consequence

The current five-point matter route has reached a clean boundary:

- the commutative algebra is forbidden at arbitrary multiplicity;
- the crossed-product algebra is selected but its canonical carriers fail;
- arbitrary proper supports exist abundantly but violate the exact
  conjugation and are not selected.

Further numerical searches inside these 256 supports would answer only
whether something can be fitted.  They cannot advance the theory unless a
new geometric principle selects a support before its Hessian behavior is
examined.
