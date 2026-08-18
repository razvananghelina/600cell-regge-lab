# Preregistration: can a second canonical axis orbit break the residual `D5`?

Date: 2026-08-11

## Disclosed starting point

Commit `a6b7e9a` proves that a six-point Hopf vacuum has stabilizer `D5` and
that this residual symmetry prevents connectedness on the fixed 936-state
carrier, even when every legal `D5` intertwiner is granted.

The same exact icosahedral geometry already contains three canonical
symmetry-axis orbits:

```text
fivefold / binary C10:   6 unoriented axes, 12 oriented points,
threefold / binary C6:  10 unoriented axes, 20 oriented points,
twofold / binary C4:    15 unoriented axes, 30 oriented points.
```

The previously certified sixth-order invariant has the fivefold orbit as
global maxima and the threefold orbit as global minima.  With one sign it can
conditionally select one orbit or the other, not both.  The twofold orbit is
only an intermediate critical orbit.  No current certified action contains
two independently fluctuating copies of this field.

The candidate observation is that a pair of axes of different orders might
have trivial common stabilizer and hence remove the residual-symmetry
obstruction.  No pair-stabilizer census has been computed before this
protocol.

## Frozen exact group construction

Construct `A5` as the 60 even permutations of five symbols.  Without using
floating coordinates, enumerate:

- all cyclic subgroups `C5`, `C3`, `C2`;
- their normalizers `D5`, `D3` and `V4`;
- the homogeneous `A5` sets `A5/C5`, `A5/C3`, `A5/C2` of oriented axis
  points;
- the homogeneous sets `A5/D5`, `A5/D3`, `A5/V4` of unoriented axes.

Verify the expected set sizes and stabilizer orders exactly.

## Frozen complete pair census

For each of the three cross-type pairs

```text
fivefold x threefold,
fivefold x twofold,
threefold x twofold,
```

and separately for oriented and unoriented versions, enumerate **every**
ordered cross-product pair.  Under the diagonal `A5` action, record before
any Dirac test:

1. total pair count;
2. complete orbit-size multiset;
3. complete common-stabilizer-order multiset, both by orbit and by point;
4. number and fraction of pairs with trivial common stabilizer;
5. number of distinct free orbits;
6. whether the action on each free orbit is a regular `A5` torsor.

No favorable relative position may be selected after this census.  Group
theory alone may distinguish double-coset orbits but does not identify a
physical incidence orbit; any later use of one orbit must cite an independent
cell-incidence relation.

## Action-provenance audit

Independently inspect the authoritative selector files and record only the
following source-level questions:

1. does one field value ever select both a fivefold and a threefold axis;
2. is a second independent vector/projector field defined;
3. is an orientation of either selected axis dynamically fixed;
4. is the relative sign that selects each orbit derived rather than assumed.

Absence is scoped to the frozen files inspected; it is not a theorem against
future dynamics.

## Acceptance and interpretation boundaries

- **DERIVED COMPLETE-BREAKING POSSIBILITY:** at least one fully enumerated
  pair orbit has trivial stabilizer.  This removes the representation-theory
  obstruction but does not select a Dirac tensor.
- **STRUCTURAL ONLY:** complete breaking requires two field values, an axis
  orientation or a relative-position orbit not selected by the existing
  action.
- **DERIVED PAIR NO-GO:** every cross-type pair retains a nontrivial common
  stabilizer.

Even a free orbit supplies a regular `A5` field.  At one vacuum, equivariance
then permits the full rectangular space `Hom_R(V_i,V_j)`, of dimension
`n_i n_j`; it does not by itself choose a matrix.  This fitting dimension and
the number of free orbits must be reported before any connectedness test.

No matter character, mass, coupling or Standard-Model target will be used.
Only a targeted verifier will be run; the full suite remains excluded by
user instruction.
