# A second unoriented axis usually does not remove the residual symmetry

Date: 2026-08-11

Protocol commit: `8384d59`.

Registered verifier:
`reproducible/verify_hopf_residual_symmetry_breaking.py`.
Targeted exact result: `14/14`.

No Dirac connectedness, matter character, mass, coupling or Standard-Model
target was evaluated.  The full suite was not run, by explicit user
instruction.

## Headline

The complete cross-orbit census refutes the simplest proposed repair.

Adding one **unoriented** threefold or twofold axis to the selected
unoriented Hopf/fivefold axis never breaks `A5` completely:

```text
D5 x D3:  60/60 pairs retain C2,
D5 x V4:  90/90 pairs retain C2.
```

Thus the fivefold Hopf vacuum plus either one of the other natural axis types
still leaves a nontrivial algebra element commuting with every covariant
Dirac operator.  The residual-symmetry connectedness obstruction survives.

There are two exact ways for a pair to have trivial stabilizer:

1. orient both axis points: every cross-type oriented pair is free;
2. abandon the fivefold axis and choose the unique free orbit among
   unoriented threefold--twofold pairs.

Neither option is selected by the current action.  This is a **DERIVED
STABILIZER CENSUS** and a **STRUCTURAL ONLY** opening, not a new finite
spectral triple.

## Exact homogeneous sets

An independent combinatorial construction starts from all 60 even
permutations of five symbols.  The exact subgroup and normalizer data are

```text
cyclic subgroup type     number     normalizer       order
C5                          6       D5                 10
C3                         10       D3                  6
C2                         15       V4                  4
```

The corresponding oriented and unoriented homogeneous sets are

```text
oriented:    A5/C5 = 12,  A5/C3 = 20,  A5/C2 = 30,
unoriented:  A5/D5 =  6,  A5/D3 = 10,  A5/V4 = 15.
```

These recover exactly the icosahedral vertex, face-centre and edge-centre
counts without floating coordinates.

## Complete unoriented pair census

| Pair | Total pairs | Orbit sizes | Common stabilizers | Free pairs |
|---|---:|---|---|---:|
| fivefold x threefold | 60 | `30,30` | `C2,C2` | `0/60` |
| fivefold x twofold | 90 | `30,30,30` | `C2,C2,C2` | `0/90` |
| threefold x twofold | 150 | `60,30,30,30` | `1,C2,C2,C2` | `60/150` |

Every orbit-size/stabilizer pair obeys exact orbit--stabilizer.  The sole
free unoriented orbit is a regular `A5` torsor of size 60.

The two fivefold--threefold orbits are the two possible relative-position
classes that remain after orientations are forgotten.  Group theory alone
does not label one as physical incidence, and neither is free anyway.

## Complete oriented pair census

All cross-type oriented pairs have trivial common stabilizer:

| Pair | Total pairs | Free `A5` orbits | Free pairs |
|---|---:|---:|---:|
| `C5 x C3` | 240 | 4 | `240/240` |
| `C5 x C2` | 360 | 6 | `360/360` |
| `C3 x C2` | 600 | 10 | `600/600` |

Each orbit has size 60 and is an exact regular `A5` torsor.  Coprime cyclic
stabilizers make trivial intersection possible, but the orbit multiplicities
remain a look-elsewhere count: incidence must select one of 4, 6 or 10
relative positions before it can carry evidence.

The old orbifold calculation does supply oriented cell incidence as one
double coset.  It does not, however, supply a dynamical orientation of the
unoriented Hopf order parameter or a lift of that incidence to the 936-state
Dirac carrier.

## The price of complete breaking

At a free vacuum the stabilizer is trivial.  Equivariance then permits

```text
Hom_1(V_i,V_j)=Hom_R(V_i,V_j),
```

whose off-diagonal dimensions are

```text
36, 72, 72, 72, 72, 144.
```

For the three legal links of a grading, the complete rectangular tensor
freedom is

```text
216 dimensions for the four edge-first readings,
180 dimensions for the four reflection-first readings.
```

This is the opposite problem from residual symmetry.  A free vacuum removes
the invariant projector obstruction, but representation theory then selects
no matrix at all.  A generic element would almost certainly connect the
carrier, yet choosing it would be precisely the fitted-Schur-data failure
that the programme forbids.

## Action provenance

The scoped authoritative audit covers

```text
verify_hopf_sixth_order_selector.py,
verify_hopf_hessian_action_origin.py,
hopf_sixth_order_selector.json,
hopf_hessian_action_origin.json.
```

It finds:

- one three-component vector field, not two independent fields;
- unoriented projectors/axis lines, with no selected sign orientation;
- one sign branch selects the six fivefold axes and the opposite sign selects
  the ten threefold axes, never both values of one field;
- the relative sign is explicitly not derived;
- the twofold orbit is intermediate, not a selected minimum.

This is a **DERIVED REPOSITORY-STATE NEGATIVE** scoped to those files.  It is
not a theorem that future dynamics cannot add a second field.

## Hostile framing audit

1. "Two different axes generically fix a frame" is not enough.  For the
   unoriented Hopf pairs, every one of 150 cases was enumerated and all retain
   `C2`.
2. Orienting an axis is physical information.  Replacing a line projector by
   a signed vector after seeing the obstruction is an extra Ising-like field,
   not a convention.
3. The one free unoriented `D3 x V4` orbit omits the Hopf/fivefold vacuum that
   motivated the six-fibration algebra.  Using it would change the vacuum
   story rather than complete it.
4. Trivial stabilizer is necessary here, not sufficient.  It removes the
   universal commutant witness but leaves 180 or 216 tensor coefficients
   unselected.
5. The incidence double coset can select a relative-position orbit, but it
   cannot by itself choose a rectangular map inside the full free-vacuum Hom
   space.

## Status ledger

- **DERIVED:** all cyclic subgroup, normalizer and homogeneous-set counts.
- **DERIVED:** complete oriented and unoriented cross-pair orbit censuses.
- **DERIVED NEGATIVE:** every unoriented Hopf--threefold and Hopf--twofold
  pair retains `C2`.
- **DERIVED:** exactly one of four unoriented `D3 x V4` orbits is free,
  containing `60/150` pairs.
- **DERIVED:** all 20 oriented cross-type orbits are free regular `A5`
  torsors (`4+6+10`).
- **DERIVED LOOK-ELSEWHERE COST:** a free vacuum exposes 180 or 216 dimensions
  of legal rectangular tensor freedom per reading.
- **DERIVED REPOSITORY-STATE NEGATIVE:** no second independent field, axis
  orientation or fixed relative sign occurs in the authoritative selector
  construction.
- **STRUCTURAL OPENING:** an oriented two-axis field or the free unoriented
  `D3 x V4` orbit could remove residual symmetry.
- **OPEN:** a geometry-selected action and incidence tensor that choose such
  a vacuum and its Dirac matrices before connectedness is inspected.
- **NO DIRAC TARGET:** connectedness was intentionally not tested in this
  preregistered enumeration stage.

## Programme boundary

The most economical completion of the existing Hopf vacuum is not merely "a
second axis".  It requires at least one additional orientation-sensitive
field, because every unoriented pair containing the Hopf axis retains `C2`.

Until such a field and its action are derived, moving to the full rectangular
free-vacuum maps would trade a clean residual-symmetry no-go for hundreds of
fittable coefficients.  That is not progress toward a selected theory.
