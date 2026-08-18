# STEP 1: blind central Krajewski census, no Hessian comparison

Date: 2026-08-11

Initial protocol commit: `8b21830`.

Weighted-minimum correction commit: `ed4b989`.

Registered verifier:
`reproducible/verify_hopf_hessian_crossed_krajewski_blind.py`.
Targeted exact result: `12/12`.

No Hessian, `Box` configuration or selector target was used.

## Arena and gates

The fixed algebra has central block sizes

```text
n=(5,5,5,15).
```

A positive-grading central support is an oriented simple graph on four nodes.
KO6 supplies the transposed negative support.  Metric-zero orientability
forbids loops and simultaneous reverse edges.  The blind binary census tests
all

```text
3^6=729
```

supports against:

1. faithful node coverage;
2. nondegenerate `Q=mu-mu^T`;
3. at least one first-order-legal odd cell block;
4. connectedness of the central commutator-link graph.

These are necessary central gates.  They do not assert a full numerical
Dirac operator or matrix-factor connectedness.

## Complete binary ledger

| gate | count |
|---|---:|
| total supports | 729 |
| faithful | 636 |
| nondegenerate Poincare form | 484 |
| nonzero first-order support | 642 |
| centrally connected | 316 |
| all four gates | 256 |

Every rank-four support is automatically faithful in this four-node census.
The 256 survivors have Hilbert-dimension multiset

| total KO6 dimension | number |
|---:|---:|
| 300 | 24 |
| 400 | 48 |
| 450 | 48 |
| 500 | 24 |
| 550 | 48 |
| 600 | 64 |

Their intersection determinants are

```text
det(Q)=1 : 144
det(Q)=4 :  96
det(Q)=9 :  16.
```

Thus abstract existence is abundant, not rare:

```text
binary hit fraction = 256/729.
```

After quotienting by all abstract permutations of the three equal `M5`
blocks and by global grading reversal, 22 support classes remain.  This is
still a substantial look-elsewhere space, and that quotient is more generous
than the actually derived symmetry because the trivial and nontrivial
stabilizer characters are distinguishable as `A5` modules.

## Correction and true weighted minimum

The initial protocol incorrectly said multiplicity could not change
intersection rank.  The correction `ed4b989` was committed before the
weighted result.

A binary survivor gives the upper bound 300.  Every integer signed
multiplicity assignment below that bound lies in the exact box

```text
[-6,6]^3 times [-2,2]^3,
13^3*5^3=274625 assignments.
```

After the dimension cut, 773 assignments remain.  Exhausting all of them
gives exactly

```text
weighted survivors <=300 = 24,
true minimum dimension   = 300,
minimum designs          = 24.
```

Every minimizer is binary.  So the original numerical minimum happens to be
correct, but only the corrected weighted enumeration proves it.

## Exact conjugation obstruction

The real geometry has the exact character conjugation swapping the two
nontrivial `M5_+` and `M5_-` nodes.  On the four central nodes this is one
transposition, represented by `P` with `det(P)=-1`.

For every alternating four-by-four intersection form,

```text
Pf(P Q P^T)=det(P) Pf(Q)=-Pf(Q).
```

Both possible grading behaviors fail:

- if conjugation preserves grading, `PQP^T=Q`, hence `Pf(Q)=0`;
- if conjugation reverses grading, `PQP^T=-Q`, but in four dimensions
  `Pf(-Q)=Pf(Q)`, again forcing `Pf(Q)=0`.

The verifier checks the identities symbolically and solves both generic
linear constraint spaces.  This is independent of multiplicity bounds.

Correspondingly, among the 256 binary survivors:

```text
conjugation-preserving grading survivors = 0
conjugation-reversing grading survivors  = 0.
```

## Blind verdict

- **STRUCTURAL EXISTENCE:** 256 abstract binary supports pass the necessary
  central gates.
- **SELECTION NEGATIVE:** there are 24 minimum designs and 22 very generously
  quotiented support classes, not a unique carrier.
- **DERIVED CONJUGATION NO-GO:** no nondegenerate intersection form with
  arbitrary multiplicities is compatible with the exact `chi<->chibar`
  conjugation, whether it preserves or reverses grading.
- **NO TARGET COMPARISON:** no independently symmetry-selected candidate
  survives, so STEP 2 is not licensed.

## Scope

The 256 supports are not full spectral triples.  They pass only the stated
central combinatorial gates.  Conversely, the conjugation Pfaffian argument
is a genuine exact obstruction under its complete symmetry hypothesis.
Abandoning that conjugation produces many abstract designs but also abandons
the only current geometric rule capable of choosing among them.
