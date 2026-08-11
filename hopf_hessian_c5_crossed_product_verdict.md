# The five-point crossed product selects an algebra, but not a viable carrier

Date: 2026-08-11

Protocol commit: `c3c1584`.

Target-blind STEP 1 commits: `c212e2d`, amended by `af05ff0` before any target
comparison.

Registered verifier:
`reproducible/verify_hopf_hessian_c5_crossed_product_blind.py`.
Targeted exact result: `17/17`.

## Headline

**STRUCTURAL ALGEBRA ADVANCE:** the five-point geometry selects the exact
noncommutative algebra

```text
C(P) crossed_product A5
  ~= M5(C)+M5(C)+M5(C)+M15(C).
```

**DERIVED NEGATIVE FOR THE CANONICAL CARRIERS:** none of the representations
or bimodules furnished without an additional Krajewski choice reaches the
selector-comparison gate.  Consequently no Hessian target comparison was
performed.

This is not a no-go for every abstract bimodule over the crossed product.  It
is a negative answer to the evidentially stronger question: does the
canonical algebra construction also select a usable physical carrier?

## Exact algebra derivation

The derived action on five primitive projectors is faithful and transitive.
A point stabilizer is the exact subgroup

```text
V4 normal in A4,       |A4|=12,
class sizes            1,3,4,4,
A4/[A4,A4]             =C3,
complex irrep degrees  1,1,1,3.
```

The transitive-groupoid theorem then gives

```text
B ~= M5(C[A4]),
Wedderburn sizes=(5,5,5,15),
dim_C B=300,
dim_C Z(B)=4.
```

This type is geometry-derived, not fitted.  It is also plainly not the
Standard-Model algebra type.

## What the physical five-dimensional modules see

The verifier constructs all three one-dimensional stabilizer characters in
exact `Q(omega)` arithmetic.  Their induced covariant representations obey
the group law and diagonal covariance exactly:

| branch | restriction to `A5` | image of `B` | kernel |
|---|---|---:|---:|
| trivial | `1+4` | full `M5` | 275 |
| `chi` | irreducible `W5` | full `M5` | 275 |
| `chibar` | irreducible `W5` | full `M5` | 275 |

Thus the crossed product does not refine the two physical monomial systems
into a smaller algebra.  It maximizes each natural image to `M5` while
discarding 275 algebra dimensions.

Every ten-state double made from two of the three branches fails blindly:

- same branch: the graph image is full `M5`, so order zero fails;
- different branches: a central projection separates the sheets, so first
  order forces every nonzero grading-odd block to vanish;
- every branch is nonfaithful to the full crossed product.

The ordered branch census is `N=9`: three order-zero failures and six
first-order failures.  The physical `W5/W5` subcensus has `N=4` with the same
two-versus-two split.  No matrix entry of the Hessian is used in this result.

## Faithful canonical carriers

A faithful left representation must contain one fundamental representation
of every simple block, so its minimum dimension is

```text
5+5+5+15=30.
```

The standard graded double has dimension 60, but its commutant has dimension
only

```text
4*(2^2)=16.
```

Order zero would embed the faithful 300-dimensional opposite algebra in that
commutant, which is impossible.

The regular bimodule has dimension 300 and only the four diagonal central
cells.  The full enveloping bimodule has dimension `300^2=90000` and all
sixteen central cells.  Their standard odd doubles represent every
Hochschild zero-cycle with the same sheet profile `(1,1)`, while the grading
has `(1,-1)`.  Both therefore fail metric-dimension-zero orientability.

These three failures cover the carriers supplied functorially by the algebra
or by its complete irreducible census.  Choosing a proper subset of the
sixteen central bimodule cells would be new Krajewski data.

## Why STEP 2 was not run

The preregistered rule allowed a Hessian comparison only after a canonical
carrier passed faithfulness, order zero and orientability.  The number of
survivors is zero.  Reporting a target hit fraction for rejected carriers
would therefore be meaningless.

This is stronger provenance than finding a convenient abstract support and
then calling it canonical.  The crossed product fixes the Wedderburn blocks;
it does not fix one of their many possible bimodule graphs.

## Status ledger

- **DERIVED:** `B=M5+M5+M5+M15`, dimension 300 and centre dimension four.
- **DERIVED:** the two nontrivial induced branches are the exact conjugate
  monomial `W5` systems.
- **DERIVED NEGATIVE:** every natural five-dimensional image is nonfaithful
  and equal to full `M5`.
- **DERIVED NEGATIVE:** zero of nine natural ten-state branches passes the
  cheap real-triple gates.
- **DERIVED NEGATIVE:** the 60-state minimum faithful left double fails order
  zero by `16<300`.
- **DERIVED NEGATIVE:** the canonical regular and full-enveloping odd doubles
  fail metric-zero orientability.
- **STRUCTURAL ADVANCE:** the noncommutative algebra type is selected by the
  geometry.
- **STRUCTURAL EXISTENCE / SELECTION NEGATIVE:** the complete blind census of
  proper central supports is recorded in
  `hopf_hessian_crossed_krajewski_verdict.md`.  It finds `256/729` necessary-
  gate survivors, 24 minimum designs, and no nondegenerate support compatible
  with exact character conjugation.
- **OPEN PHYSICS:** a geometry-selected faithful real bimodule and a licensed
  coupling of the Hessian selector to it.

## Programme consequence

The two canonical five-point algebra choices are now bounded:

1. `C^5`: no nonzero `A5`-equivariant KO6 metric-zero-orientable bimodule
   exists, even with arbitrary multiplicities.
2. `C(P) crossed_product A5`: the algebra is noncommutative and canonical,
   but none of its canonical carriers is a viable finite triple.

The abstract central-cell enumeration has now been performed blindly.  It
answers existence positively but selection negatively.  No independent
geometric rule chooses among the 256 survivors; the exact available
conjugation instead obstructs all of them at Poincare rank.  Consequently no
Hessian comparison is licensed in this arena.
