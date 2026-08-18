# Orbifold groupoid and KO6 Poincare-duality audit

Date: 2026-07-28

## Decision

Two decisive negative results were obtained.

1. The full-character incidence construction on the canonical
   `C10/C4/C6` orbifold groupoid has intersection-form rank `4/20`.
   Poincare duality fails.  **DERIVED negative for the specified canonical
   construction.**
2. In KO-dimension six the intersection form is antisymmetric.  Therefore
   any candidate algebra whose relevant `K0` has odd rank is intrinsically
   Poincare-degenerate.  In particular the nine-block full algebra
   `C[2I]` is excluded before choosing `D`.  **DERIVED no-go.**

No Dirac operator or grading was fitted.

## 1. The KO6 parity lemma

Let `p_i` be self-adjoint projections representing a basis of `K0(A)` and
put

`q_i=J pi(p_i) J^-1`.

Assume order zero, `[gamma,pi(A)]=0`, `J^2=1`, and the KO6 relation

`J gamma J^-1=-gamma`.

Define

`cap_ij=Tr(gamma pi(p_i) q_j)`.

Conjugation by the antiunitary `J` gives

`conj(cap_ij)=Tr(J gamma pi(p_i)q_j J^-1)`

`             =-Tr(gamma q_i pi(p_j))`

`             =-cap_ji`,

where order zero permits interchange of `q_i` and `pi(p_j)`.  The trace is
real because `gamma`, `pi(p_i)`, and `q_j` are commuting self-adjoint
operators.  Hence

`cap_ji=-cap_ij`.

The intersection form is antisymmetric.  An antisymmetric matrix of odd
order has zero determinant.  **DERIVED.**

### Consequences

- `C[2I]` has nine complex simple blocks, so its standard `K0` rank is nine:
  Poincare duality is impossible in KO6.  **DERIVED negative.**
- The literal three-summand `C+H+M3(C)` corner also has odd `K0` rank and
  cannot by itself carry a nondegenerate KO6 pairing under this axiom set.
  It was already nonunital on the arena.  **DERIVED negative, scoped to the
  literal three-summand algebra.**
- The orbifold groupoid below has 20 simple blocks, so it passes this parity
  screen.  Even rank is necessary, not sufficient.

This lemma does not claim that every presentation called a Standard-Model
finite geometry uses precisely this unmodified algebra and Poincare axiom.
Extra summands, reductions, quotients, or weakened axioms must be audited
separately.  **OPEN outside the stated hypotheses.**

## 2. Canonical transformation groupoid

The three derived transitive `2I` sets are

`2I/C10`, `2I/C4`, `2I/C6`,

of sizes `12,30,20`.  For a transitive finite action,

`C[G lt G/H] = M_[G:H](C[H])`

up to canonical groupoid-algebra isomorphism.  Since all stabilizers are
cyclic, the three components decompose as

`10 M12(C) direct-sum 4 M30(C) direct-sum 6 M20(C)`.

The dimension is

`10*12^2 + 4*30^2 + 6*20^2 = 7440`.

It has 20 simple blocks and is Morita equivalent to

`C[C10] direct-sum C[C4] direct-sum C[C6] = C^20`.

All statements in this section are **DERIVED** from the certified orbit and
stabilizer data.

## 3. Canonical incidence correspondence

There are 60 vertex-edge flags and 60 edge-face flags.  `A5` acts
transitively on each flag set, so its flag stabilizer is trivial.  Its
inverse image in `2I` is the central `C2`.

Every cyclic-stabilizer character restricts to one of the two characters of
this `C2`.  Frobenius reciprocity gives a one-dimensional flag intertwiner
exactly when the two exponent parities agree.  Therefore:

- the `C10-C4` compatibility matrix is `10 x 4`, has 20 nonzero entries,
  and rank 2;
- the `C4-C6` matrix is `4 x 6`, has 12 nonzero entries, and rank 2;
- the canonical KO6 antisymmetrization on the 20 `K0` generators has rank
  exactly 4 and nullity 16.

Uniform nonzero weights on the two flag orbits do not change the rank.
**DERIVED.**

There is a stronger boundary.  Cellular parity gives 16 even generators
(10 vertex plus 6 face) and only 4 odd edge generators.  Any pairing built
only from adjacent-degree correspondences has block form

`[[0,B],[-B^T,0]]`, with `B` of size `16 x 4`.

Its rank is at most 8, even if arbitrary character-dependent multiplicities
are allowed.  Thus fitting weights cannot repair Poincare duality inside
this adjacent cellular design.  **DERIVED no-go under the cellular
incidence boundary.**

## 4. What is structural

Using every stabilizer character is the canonical regular twisted-sector
completion once the transformation groupoid is selected.  Selecting the
transformation groupoid itself as the internal algebra remains a
**STRUCTURAL** theoretical move; it is not forced by a previous theorem.

The Krajewski interpretation of cellular incidence uses the derived
cell-degree grading and flag correspondences.  It is the minimal canonical
candidate, not a classification of every possible bimodule over the
20-block groupoid algebra.  **STRUCTURAL candidate with a DERIVED negative
verdict.**

## 5. Status ledger

### DERIVED

- KO6 antisymmetry of the intersection form;
- odd `K0` rank implies Poincare degeneracy;
- exclusion of the nine-block full `C[2I]` algebra;
- exact groupoid Wedderburn decomposition and dimension 7440;
- central `C2` flag stabilizers and parity compatibility;
- canonical incidence pairing rank `4/20`;
- rank ceiling `8/20` for all adjacent-degree cellular multiplicities.

### STRUCTURAL

- promoting the transformation groupoid to the internal algebra;
- taking the complete regular character fibers as the twisted-sector
  Hilbert input;
- interpreting the cellular flag correspondence as the matter Krajewski
  diagram.

### PATTERN

- none is promoted to particle physics;
- the appearance of 20 `K0` sectors is not identified with a multiplet
  count.

### OPEN

- a geometrically forced even-rank quotient/subalgebra with a nondegenerate
  pairing;
- a non-cellular correspondence using more than adjacent flag data;
- an inter-level Bratteli/Kasparov class whose pairing is nondegenerate;
- any Standard-Model algebra, hypercharge, generations, or Yukawa sector.

## 6. New surviving structure: the oriented chamber double

The barycentric subdivision contains complete flags

`vertex < edge < face`.

There are

`20*3*2=120`

such chambers.  The orientation-preserving rotation group `A5` has order
60, so the chambers split into exactly two free `A5` orbits.  They are the
two orientations of a barycentric chamber.  Their lifts have central `C2`
stabilizer in `2I`.

The quotient orbifold cell count closes exactly:

`1/10+1/4+1/6 - 3/2 + 2/2 = 1/60 = chi(S2)/120`.

Thus the two chamber sheets are not an inserted copy.  They are a
**DERIVED geometric double**.  Orientation reversal exchanges them.

It does **not** implement the outer/Galois automorphism.  The full geometric
icosahedral symmetry group is `A5 x C2`; an improper symmetry is central
inversion times a rotation and induces only an inner automorphism on `A5`.
The proposed geometric identification with the Galois outer involution is
therefore **REFUTED**.

This does not repair the rank-4 pairing computed above, because that pairing
used only the minimal adjacent stabilizer-character correspondence.  It
does identify the next nontrivial experiment: construct the full
barycentric chamber representation, retain its two orientation orbits, and
test whether its Galois/real pairing avoids sheet cancellation.

## Final boundary

The stabilizer-groupoid idea does not open the matter gate in its canonical
minimal form.  Its failure is informative: local flag transport collapses
all stabilizer characters to the two central-parity channels.  A successful
construction must contain genuinely richer information than stabilizer
restriction along one icosahedral cell incidence.
