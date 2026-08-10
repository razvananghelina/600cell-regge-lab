# Chamber rigidity after the noncommutative counterexample

Date: 2026-08-09

> **2026-08-10 resolution:** B1 is now **DERIVED REFUTED** by an exact
> `M2(C)+C^3` representation satisfying orientability, connectedness, and all
> other stated gates on the fixed chamber matrices. See
> `chamber_b1_refutation.md`. This note retains the earlier A2 audit and the
> historical proof routes that preceded that witness.

## Complete statements under audit

### A2: partition-amplification conjecture

Fix the derived oriented-chamber carrier `H=C^120`, adjacency `D`, grading
`gamma`, and reflection real structure `J`.  Let `A_P` be a faithful unital
partition algebra satisfying first order, metric-dimension-zero
orientability, spectral connectedness, and nondegenerate integral
intersection form.  The conjecture asserted that the graph on central blocks

`i ~ j iff c_ij=dim(P_i J(P_j)H_+)=1`

must be connected.

**Verdict: REFUTED.**

### B1: repaired unrestricted commutativity conjecture

On the same fixed `H,D,gamma,J`, suppose a finite-dimensional unital star
algebra has a faithful representation satisfying:

1. `[gamma,A]=0`;
2. order zero;
3. first order;
4. nonzero represented one-forms;
5. a nondegenerate intersection form;
6. metric-dimension-zero orientability;
7. connectedness, `[D,a]=0` only for scalar `a`.

Does this force `A` to be commutative?

**Verdict: DERIVED REFUTED (2026-08-10).** The later pure-source witness
constructs a noncommutative algebra satisfying all seven gates for the fixed
`D`; see `chamber_b1_refutation.md`.

## A1: census of the witnesses already in hand

The fixed-`C5` contraction action has 18 affine free bits.  Reenumerating all
`2^18=262144` vectors reproduces the declared set of 84 candidates whose
determinants are then individually recomputed exactly over `Z`.

The unit-cell connectivity distribution is:

- 72 connected;
- 12 disconnected.

Every witness has algebra dimension 30.  The complete component-size
distribution is

- 72 with `(30)`;
- 4 with `(6,6,6,6,6)`;
- 8 with `(4,4,4,4,4,2,2,2,2,2)`.

The number of unit cells has distribution

`{40:8, 50:4, 60:44, 70:4, 80:24}`.

Each retained witness has an exact determinant certificate.  The discovery
stage still uses numerical rank as a prefilter, so 84 is the exact size of
the declared discovered set, not a certified statement that no additional
unimodular vector was missed by a numerical false negative.  **DERIVED for
the declared 84; exhaustive-total count not claimed.**

For the older `C^36` certificates:

- PALPABLE: `det=4`, 100 unit cells, connected unit graph;
- INTEGRAL: `det=1`, 104 unit cells, connected unit graph.

PALPABLE is full-rank but not strict integral PD.  It is included because A1
explicitly requested both registered certificates.

## Exact A2 counterexample

One disconnected witness is encoded by

```text
011001110110100111001000101001100011001101101011010000110101101110000111100101110010011001
```

It has:

- `A=C^30`;
- first order on every chamber edge;
- metric-dimension-zero orientability;
- connected quotient graph, hence spectral connectedness;
- `det(Cap)=1` and an integral Pfaffian certificate;
- unit-cell components `(6,6,6,6,6)`.

Thus all hypotheses of A2 hold and its conclusion fails.  **DERIVED exact
refutation.**

## The corrected amplification lemma

Connectivity of the unit-cell graph was stronger than necessary.  If a node
`i` is incident to any unit cell `(i,j)`, order zero for a proposed matrix
amplification requires

`n_i*n_j divides c_ij=1`.

Hence `n_i=n_j=1`.  It is sufficient that the unit-cell graph have no
isolated vertices; its distinct components need not communicate.

All 84 `C5` witnesses and both `C^36` certificates have no isolated unit-cell
node.  Therefore every matrix amplification retaining any one of these fixed
central-support systems is scalar on every node.  **DERIVED for this finite
enumerated collection.**

This is not a theorem for all legal partition algebras.  Unimodularity alone
only makes the gcd of each intersection-form row equal to one; a row may, in
principle, contain entries such as 2 and 3 without containing a unit entry.
Whether chamber combinatorics forbids an isolated unit-cell node in every
legal partition remains **OPEN**.

## B1: locality reduction is false

Let

`S=(D J)|H_+`.

Its exact integer matrix has a six-dimensional `-1` eigenspace.  Select an
exact rational eigenvector crossing more than one block of the connected
`C5` partition and form the Householder operator

`U=I-2 vv^T/(v^T v)`.

Then `U` is rational, orthogonal, non-monomial, and commutes with `S`.  In the
`J`-paired sheet basis, `diag(U,U)` preserves `D,J,gamma`.  Conjugating the
exact `C5` partition triple by this unitary preserves order zero, first order,
orientability, connectedness, one-forms and PD, while a central indicator
projector acquires nonzero off-diagonal chamber entries.

Therefore orientability plus connectedness do not force central projectors
to be local chamber indicators.  The proposed route reducing B1 to A2 is
**DERIVED REFUTED**. The conjugated algebra remains commutative, so this does
not by itself refute B1. B1 is refuted separately by the later exact witness.

## B1: noncommutative design survivor

The cheap multiplicity filters also do not prove B1.  For

`A=M2(C) direct-sum C^3`, `n=(2,1,1,1)`,

take the signed upper-triangular multiplicities

`(-5,10,-1,1,-9,18)`

in pair order `(01,02,03,12,13,23)`, where a negative value reverses the
ordered cell.  This gives positive-sheet cells

```text
(1,0) x5, (0,2) x10, (3,0) x1,
(1,2) x1, (3,1) x9, (2,3) x18.
```

Their weighted dimension is 60.  There are no diagonal or reverse-paired
cells, so the multiplicity design is orientable.  Its antisymmetric
intersection matrix has `Pf=-1`, `det=1`. The exact capacitated matching of
Krajewski-allowed cell supports has structural rank 60. This is only the
support part of first order: a non-scalar shared factor also requires an
identity intertwiner.

This is only a **STRUCTURAL design-filter survivor**.  Structural rank is a
necessary generic invertibility condition, not a certificate that an allowed
Dirac lies in the unitary-congruence/Takagi class of the exact chamber `S`.
No fitted `D` is registered.

An alternating-projection detector was tried on this inverse-spectral gate.
Its best completed relative singular-value residual was approximately
`5.06e-3`; this is **PATTERN/INCONCLUSIVE**, neither an existence witness nor
a no-go.  Exact fixed-`D` compatibility remains **OPEN**.

## Status ledger and evasion boundaries

- **DERIVED NEGATIVE:** A2 is false; unit-cell connectivity is not forced by
  the stated partition gates.
- **DERIVED POSITIVE:** absence of isolated unit-cell nodes forces scalar
  amplification; it closes all 84 declared `C5` witnesses and both registered
  `C^36` certificates at fixed central supports.
- **DERIVED NEGATIVE:** orientability and connectedness do not force chamber-
  local central projectors.
- **STRUCTURAL:** the noncommutative `M2+C^3` multiplicity survivor.
- **PATTERN/INCONCLUSIVE:** the alternating-projection residual.
- **DERIVED NEGATIVE:** B1 is false; see `chamber_b1_refutation.md`.
- **OPEN:** the strict unimodular strengthening of B1; a theorem excluding
  isolated unit nodes for every legal partition.
- **NOT CLAIMED HERE:** the historical design survivor is not itself a full-
  gate witness. The later exact construction is one, but no Standard Model
  algebra or geometry-selected physical gauge sector is claimed.

The later B1 witness is the first certified noncommutative all-gate chamber
algebra, so the old “all witnesses are tori” statement is no longer true. It
is structurally chosen and has no derived color/matter sector, so it still
does not open the Standard-Model physics gate.

Exact verifier: `reproducible/verify_chamber_rigidity_audit.py`.
