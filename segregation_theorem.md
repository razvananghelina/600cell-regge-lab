# Central-parity segregation theorem and its exact escapes

Date: 2026-07-23

## Decision

**DERIVED:** the common roof over the four equivariant matter failures is the
central element `z=-1` of the binary icosahedral group `2I`.  If `gamma` is
its integer-spin/spinor grading, every `2I`-equivariant linear operator
commutes with `gamma`.  Consequently a `2I`-equivariant gamma-odd operator is
zero.

**DERIVED:** `2I` has exactly one involution, `-1`.  A subgroup omits `-1`
exactly when it has odd order, and the complete list of odd-order subgroup
types is

`C1, C3, C5`.

There are ten actual `C3` subgroups and six actual `C5` subgroups.  Therefore
the maximal possible residual symmetry of a nonzero parity-odd Dirac sector
is `C3` or `C5`.

On the multiplicity-one node module `W`, the exact complex dimensions of the
full parity-odd operator spaces are

| residual group | `dim_C OddEnd_H(W)` |
|---|---:|
| `C5` | 88 |
| `C3` | 148 |
| `C1` | 448 |

The same numbers are the real dimensions of the self-adjoint odd parts.

**DERIVED further combinatorial test:** after restricting the canonical
`C + H + M3(C)` seeds, `C3` admits an explicit `16+14` Krajewski reblocking
and a nonzero residual-`C3`-equivariant odd block type allowed by the
Krajewski shared-index rule.  `C5` admits no grading-wise
compatible reblocking at all.  Thus the smallest node test leaves `C3`, but
not `C5`, alive.

This is not yet existence of a real finite spectral triple or of a Dirac
meeting a matrix-level double-commutator test.  The verifier constructs the
restricted bimodule character and checks the Krajewski shared-index rule plus
the existence of a common `C3` character.  It does not construct a faithful
real-algebra representation, opposite action, `J`, or a matrix `D`.

This does **not** select a Dirac, hypercharge, or generations.  Symmetry
breaking enlarges the odd space from zero to dimension 148 in the surviving
`C3` case.  Selecting a point in that space remains **OPEN**.

All finite claims are checked exactly by
`reproducible/verify_segregation_theorem.py`.

## 1. Central-character lemma

Let `V` be a finite-dimensional real or complex representation of `2I`, let
`z=-1`, and put `Gamma=rho(z)`.  Since `z` is central and `z^2=1`,
`Gamma^2=1`; in characteristic zero,

`V = V_+ direct-sum V_-`, where `Gamma|V_+ = +1` and
`Gamma|V_- = -1`.

On an irreducible real or complex summand, Schur's lemma says that the central
element acts in the center of the relevant real division algebra.  The
additional equation `Gamma^2=1` forces the scalar to be `+1` or `-1`.
Equivalently, the real proof needs no algebraic closure: the polynomial
`x^2-1` has distinct real roots, so the two eigenspaces above are real and
`2I`-stable.

If `T` is `2I`-equivariant, then

`T Gamma = T rho(z) = rho(z) T = Gamma T`.

Hence `T(V_+)` lies in `V_+` and `T(V_-)` lies in `V_-`.  If also
`T Gamma = -Gamma T`, then `2 Gamma T=0`, and therefore `T=0`.

This proof applies to real-linear or complex-linear `T`.  It also applies to
complex anti-linear equivariant maps because the two central scalars are
real.  No assumption that all irreducibles are of real type is being made;
quaternionic spinor irreducibles still have central scalar `-1`.

**Segregation theorem (DERIVED).**  Every equivariant operator preserves the
central-character sectors.  Therefore no nonzero parity-odd Dirac can retain
full `2I` equivariance.

## 2. One theorem, four scoped corollaries

The lemma identifies the same zero intersection in the four previous
structures.  It does not erase their stronger, independent algebraic
obstructions.

1. **Node space.**  On `W` the grading is `rho_i(-1)`, with dimensions
   `16+14`.  Thus every `D in End_2I(W)` is block diagonal and every
   equivariant odd `D` is zero.  On `W+W^sigma`, the same statement holds
   isotypic block by isotypic block.  The earlier missing-fusion-support and
   `M3` obstructions remain independently valid.

2. **Bidirected edge/Hom space.**  Use the source copy of `2I` on
   `Hom(V_s,V_t)`.  Its central element has sign `rho_s(-1)`, which is exactly
   the source-parity/orientation grading used in the edge note.  A source-
   equivariant operator preserves that sign, so an equivariant odd `D`
   vanishes.  The endpoint first-order theorem was stronger: it killed all
   odd blocks for the full endpoint algebra, whether or not an arbitrary
   block was equivariant.

3. **Preprojective regular fiber.**  On `C[2I]`, gamma is left multiplication
   by `-1`.  Every left-equivariant/right-convolution operator commutes with
   gamma, reproducing the zero odd part of the existing operator family.
   The maximal diagonal-bimodule first-order obstruction and the canonical
   `J gamma=+gamma J` sign remain additional negatives.

4. **Rooted Bratteli/McKay tower.**  Tensor floor `n` has central character
   `(-1)^n`.  Any `2I`-equivariant map between consecutive floors would
   connect opposite central characters and is therefore zero.  The earlier
   cumulative-endpoint calculation additionally showed that a non-equivariant
   nonzero shift fails its chosen first-order condition.

Thus the lemma exactly reproduces `D=0` for the equivariant odd sector of all
four constructions.  Claims about their larger non-equivariant candidate
spaces continue to rest on the separately verified first-order or support
arguments.  **DERIVED scope boundary.**

## 3. Exact subgroup classification

Use the concrete group structure

`2I isomorphic to SL(2,5)`.

Enumeration of the determinant-one `2 by 2` matrices over `F5` gives 120
elements.  Solving `g^2=1` exactly gives only `1` and `-I`; hence `-I` is the
unique involution.  This is a group-structure verification, not an appeal to
the general binary-polyhedral statement.

If a finite subgroup `H` has even order, Cauchy's theorem gives an
involution, necessarily `-1`.  Conversely, if `-1 in H`, then `H` has even
order.  Therefore

`-1 notin H iff |H| is odd`.

An odd-order subgroup injects through
`2I -> 2I/{+/-1} isomorphic to A5`.  Its order divides both 120 and 60, so
the only possibilities are `1,3,5,15`.  Groups of prime order are cyclic.
Order 15 is impossible: its Sylow `C5` is unique and hence normal, so the
whole subgroup would lie in its normalizer.  Exact enumeration gives a
normalizer of order 20 in `2I`, hence order 10 in `A5`, too small to contain
an order-15 subgroup.  The enumeration also gives ten `C3` and six `C5`
subgroups.

**Escape-classification theorem (DERIVED).**  A nonzero parity-odd Dirac can
be equivariant only under `C1`, `C3`, or `C5`.  The maximal residual groups
are cyclic `C3` and `C5`.

This is a necessary bound.  It does not say that every geometric structure
admits a Dirac with either maximal residual group.

## 4. Exact first consequences on `W`

Write `W=W_+ direct-sum W_-`, with dimensions `16` and `14`.  Restriction of
the character table gives cyclic weight multiplicities

| group | `W_+` weights | `W_-` weights |
|---|---|---|
| `C3` | `(6,5,5)` | `(4,5,5)` |
| `C5` | `(4,3,3,3,3)` | `(2,3,3,3,3)` |
| `C1` | `(16)` | `(14)` |

Consequently

`dim_C OddEnd_H(W) = 2 sum_chi m_+(chi)m_-(chi)`,

which gives `148`, `88`, and `448`, respectively.  For self-adjoint odd
operators the lower-left block is the adjoint of the upper-right block, so
the real dimension is the same number.

### Canonical seed and first-order test

Keep the five Krajewski labels
`(1,1bar,2,3,3bar)` and restrict the derived scalar, weak, and color seeds to
the residual cyclic group.  Exhaustive exact dynamic programming over all 25
cells gives:

- **`C5`: DERIVED negative.**  Neither `W_+` nor `W_-` has any compatible
  nonnegative Krajewski reblocking.  The new obstruction is exact cyclic
  character content, not lack of odd equivariant operators: the latter space
  has dimension 88.
- **`C3`: DERIVED positive combinatorial candidate.**  One explicit witness is

  `H+ = (2,2) + 2(3bar,2)`,

  `H- = (2,1bar) + 2(3bar,2)`.

  Their `C3` weights are exactly `(6,5,5)` and `(4,5,5)`.  The block type
  `(2,2) <-> (2,1bar)` shares its left weak label, is first-order legal, and
  is odd in the Krajewski combinatorial sense.  The two restricted cell
  representations share nonzero `C3` character support, so a nonzero
  residual-equivariant intertwiner exists.  The witness uses scalar,
  quaternionic weak, and color seed labels.  Matrix-level first order for a
  specified real representation and opposite action remains open.

For this witness the gamma-preserving algebra commutant before choosing `D`
is

`C + M2(C) + C + M2(C)`,

one scalar and one doubled-block multiplicity factor in each chirality.
Its self-adjoint part has real dimension `1+4+1+4=10`.  It therefore admits
many possible commuting `Y` candidates, but neither the residual symmetry nor
first order selects one.  A chosen nonzero `D` can reduce this commutant in a
choice-dependent way.

**OPEN:** select `D` inside the 148-real-dimensional self-adjoint odd
`C3`-equivariant space; derive a preferred commutant generator `Y`; derive
color orientation and a generation decomposition.

**Route C does not yet apply.**  The witness is a single `30=16+14` module,
not a derived `M15` or `M16` generation repeated three times.  Without a
selected generation-blind `Y`, the anomaly equations have no licensed input.
The older conditional Route-C factorization remains true, but it is not
activated by this existence result.

## 5. Pattern flag

`C3` matches `N_gen=3`, and `C5` matches the founding integer `a1=5`.

**PATTERN only.**  The subgroup theorem derives the two maximal odd-order
possibilities, but nothing presently derives which subgroup is realized or
why its order should be identified with generations or `a1`.  Moreover, the
canonical node first-order test distinguishes them in favor of `C3`; this is
evidence about that specific structure, not yet a physical selection
principle.

## Status ledger

### Strengthened

- **DERIVED:** one central-character theorem covers the equivariant odd
  sectors of nodes, edges, the preprojective fiber, and Bratteli floors.
- **DERIVED:** `-1` is the unique involution of `2I`.
- **DERIVED:** all odd-order subgroup types are exactly `C1,C3,C5`.
- **DERIVED:** the maximal residual symmetry of a nonzero odd Dirac is `C3`
  or `C5`.
- **DERIVED:** exact odd-operator dimensions on `W` are `88,148,448`.
- **DERIVED:** the canonical `C3` seed character admits a nonzero
  residual-equivariant odd Krajewski-legal block type.
- **DERIVED negative:** the analogous `C5` grading-wise seed reblocking has
  zero solutions.

### Downgraded or delimited

- The segregation lemma covers equivariant odd operators; the older
  non-equivariant first-order no-gos remain separate stronger statements.
- Surviving `C3` combinatorial existence is not a constructed finite spectral
  triple and does not certify the first-order double commutator.
- `C3=3 generations` and `C5=a1` are **PATTERN**, not derivations.
- Route C is not triggered by the 30-dimensional witness.

### Open

- a geometric choice between the residual subgroup possibilities;
- a selection principle for `D`, quantified by the 148-dimensional `C3`
  self-adjoint odd space;
- a canonical `Y`, color orientation, and generation census;
- a physical explanation, if any, of the `3/5` pattern.
