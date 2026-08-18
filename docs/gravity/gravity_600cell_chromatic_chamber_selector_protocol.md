# Preregistration: static chamber selector for the chromatic `Z2`

Date: 2026-08-17

Prior-art gate: `e946b0a`.

Status: frozen before the sign-action enumeration.  The expected
`STATIC_CHIRAL_SELECTOR_NO_GO` outcome is disclosed in the prior-art note.

## 1. Frozen inputs and exclusions

Use only:

```text
commons/cell600.py
SHA-256 ea5bce4b6c52e0834539ca4b1df9c6a67a3a5ed4da32f4e0298a493fc5315c7f

reproducible/gravity_600cell_chromatic_cover_orbits.json
SHA-256 682e3cfaa0c2912085c0375281817e217f19a54bfc9d6ec9b296844063be7121
```

The artifact is a frozen derived input, not a result to be re-fitted.  Require
its `22/22`, `CHIRAL_COVER_AMBIGUITY`, left/right orbit, improper-exchange and
full-`H4` `Z2` fields before using them.

Independently rebuild the chamber geometry.  Do not parse the earlier chamber
or Hopf-axis result artifacts.

Exclude every Regge action, nonlinear schedule value, Standard-Model
chirality, matter character, preferred phase, mass, coupling and continuum
target.  No base chamber or hand-picked bijection between the two 120-state
sets is allowed.

## 2. Chamber reconstruction

Reconstruct the regular icosahedron and require f-vector `(12,30,20)`.  Build
all complete flags `vertex < edge < face` and require 120 chambers.

Generate the 60 proper rotations as exact vertex permutations.  Require two
free chamber orbits of size 60.  Build central inversion `J`, require that it
is a fixed-point-free involution, commutes with every rotation and exchanges
the two sheets.  Define `gamma=+1/-1` on the two sheets and require

```text
gamma(Jc) = -gamma(c).
```

The arbitrary global naming of the sheets must not enter any count.

## 3. Split versus binary 120-state group laws

Construct the full chamber permutation group from the 60 rotations and their
products with `J`.  Require a free transitive action of order 120 and compute,
from its permutation multiplication table:

- the number of nonidentity involutions;
- the order of its commutator subgroup.

Independently reconstruct the `2I` quaternion multiplication table on the
600-cell vertices and compute the same two invariants.  Unequal invariants
prove that the two regular group laws are not isomorphic.  This statement is
scoped: it forbids a group-law-preserving canonical identification, not every
possible correspondence after adding extra data.

## 4. Exhaustive two-sign map census

Use sign values `{-1,+1}`.  Under a proper symmetry all signs are fixed.  Under
reflection impose the already derived laws

```text
gamma -> -gamma,
s     -> -s,
d     -> -d,
chi=s*d -> chi.
```

Enumerate all four functions between each pair of two-point sets.  Report
both total equivariant maps and equivariant bijections for

```text
gamma -> s,
gamma -> d,
gamma -> chi.
```

Also enumerate the reflection orbits on the four-state products
`(gamma,s)`, `(gamma,d)` and `(gamma,chi)`, and the fixed-point orbits of
`chi` itself.  An invariant potential has one unconstrained value per orbit;
two orbits are nonselection, not a unique minimum.

## 5. Mechanical outcome

- `UNIQUE_STATIC_CHIRAL_SELECTOR` only if the exact symmetry data leave one
  and only one invariant value of `chi`, or one uniquely defined equivariant
  construction fixes it.
- `STATIC_CHIRAL_SELECTOR_NO_GO` if all controls pass, `gamma` is odd while
  `chi` is even, `gamma->s` and `gamma->d` each have two equivariant
  bijections, `gamma->chi` has none, and `chi` has two fixed singleton
  orbits.
- `OPEN_CONTROL_FAILURE` otherwise.

## 6. Interpretation boundary

This protocol asks only whether **existing static symmetry** selects the
residual chirality.  It does not test a dynamical parity-breaking term.  Such
a term remains admissible only after its exact form and coefficient sign are
derived without looking at which schedule behaves favourably.

Register the verifier before its first execution and run only that verifier.
The full suite remains excluded by user instruction.
