# Preregistration: does the binary/chamber geometry orient the six Hopf axes?

Date: 2026-08-11

## Disclosed motivation

The complete pair census in commit `4ad27bc` shows that every unoriented pair
containing a Hopf/fivefold axis retains a common `C2` stabilizer.  To remove
the residual-symmetry obstruction while preserving the Hopf vacuum, at least
one orientation-sensitive datum is required.

The repository already contains three objects that could be mistaken for
such a datum:

1. the two handed fibration families `qH` and `Hq`;
2. the two oriented chamber orbits;
3. the incidence-selected fibre-edge pair `{r,r^-1}` inside the `D5`
   stabilizer.

The candidate negative was noticed before implementation: an oriented
fivefold axis is the homogeneous cover

```text
A5/C5 -> A5/D5,
```

and the nontrivial element of `D5/C5` reverses it.  Thus a canonical
orientation would require an `A5`-equivariant section that may not exist.
This expectation is disclosed; the calculation is not blind.

## Frozen exact construction

Construct `A5` as all 60 even permutations of five symbols.  Enumerate all
six `C5` subgroups, their order-ten normalizers `D5`, the six-point axis set
`X=A5/D5`, and the twelve-point oriented cover `X~=A5/C5`.

Construct the quotient map `pi:X~->X` explicitly from cosets and verify:

- every fibre has exactly two points;
- the map is `A5`-equivariant;
- a normalizer reflection fixes the underlying axis, swaps the two cover
  points, and conjugates every order-five generator to its inverse class.

## Frozen section and carrier audit

Compute exactly the number of `A5`-equivariant maps from the six-axis orbit
to each already present orientation candidate:

```text
A5/C5                         oriented fivefold points,
(A5/D5) x {qH,Hq}            handed fibration labels,
(A5/1) x {chamber +,-}       oriented chamber sheets.
```

For the first and third targets, record whether any map is a section of the
natural projection back to the six axes.  For the handed double, record
whether its two equivariant choices distinguish the two points over an axis
or merely select a handed copy of the same unoriented label.

This is a homogeneous-set calculation: for `G/H -> Y`, equivariant maps are
in bijection with points of `Y` fixed by `H`.  The verifier must also exhaust
maps directly, so the conclusion does not rest on quoting the theorem.

## Frozen global-orientation census

A global fibre orientation chooses one point in each of the six two-point
fibres of `pi`, hence there are exactly `2^6=64` assignments.  Enumerate all
64 under the induced `A5` action and record:

1. the complete orbit-size multiset;
2. the complete stabilizer-order multiset, by orbit and by assignment;
3. the number fixed by all of `A5`;
4. the number with trivial stabilizer;
5. the number of free `A5` orbits;
6. whether choosing only one oriented Hopf axis leaves stabilizer `C5`.

No favorable assignment may be chosen after this census.  If free
assignments exist, their count is a look-elsewhere denominator, not a
selection theorem.

## Acceptance and kill boundaries

- **DERIVED STATIC ORIENTATION:** one of the already defined handed/chamber
  structures supplies an `A5`-equivariant section of `A5/C5 -> A5/D5`.
- **STRUCTURAL OPENING:** no invariant section exists, but some of the 64
  global assignments have trivial stabilizer.  Then an additional
  orientation field could remove the obstruction, but its vacuum remains
  unselected.
- **DERIVED STATIC-ORIENTATION NO-GO:** no existing candidate supplies a
  section and every global assignment retains a nontrivial stabilizer.

The claim is scoped to static `A5`-natural data.  A dynamical pseudoscalar,
boundary condition or chosen chamber could orient the fibres, but would be a
new symmetry-breaking field.

No Dirac matrix, matter character, mass, coupling or Standard-Model target
will be inspected.  Only a targeted verifier will be run; the full suite
remains excluded by user instruction.
