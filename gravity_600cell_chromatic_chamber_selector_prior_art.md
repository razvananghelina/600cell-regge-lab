# Prior-art gate: can chamber `gamma,J` select chromatic chirality?

Date: 2026-08-17

Status: written before the new sign-action census.  The expected no-go is
disclosed from already certified repository facts.

## Exact question and complete hypotheses

The previous all-cover audit derived three signs:

```text
s = +1/-1   for left/right 2T-coset cover chirality;
d = sign of the chromatic degree;
chi = s*d.
```

Every proper `H4` rotation preserves `s` and `d`; every improper quaternion-
conjugating symmetry reverses both.  Hence `chi` is even under reflection.

Independently, the fixed 120-state icosahedral chamber carrier has grading
`gamma=+1/-1` on its two free `A5` chamber sheets and geometric reflection
`J`, with

```text
J gamma J^-1 = -gamma.
```

The question is narrowly:

> Do the already derived static data `gamma` and `J`, without a new action
> coefficient, base chamber, section, particle target or chosen bijection,
> select one value of the residual chromatic invariant `chi`?

This is not the broader question whether a later parity-violating dynamical
term could select it.

## Relevant prior art

- Fisk proves the ten left/right-coset five-colourings of the 600-cell in
  [Coloring the 600 Cell](https://arxiv.org/abs/0802.2533).
- The binary icosahedral group is the non-split perfect central extension
  `2.A5 ~= SL(2,5)`, whereas full icosahedral chamber symmetry is the split
  group `A5 x C2`.  The non-split Schur-cover description is standard and is
  used explicitly, for example, in [Frobenius groups and retract
  rationality](https://doi.org/10.1016/j.aim.2013.06.008).
- For transitive `G`-sets, equivariant maps are determined by the image of a
  base point subject to its stabilizer.  This is the standard homogeneous-
  space fixed-point criterion.

The closest internal prior result is stronger than a mere analogy:
`hopf_axis_orientation_verdict.md` already proves that chamber orientation
does not select handed Hopf-axis data and that left/right handed copies remain
two equally natural equivariant choices.  The present mission applies the
same symmetry discipline to the newly derived `2T`-cover/degree invariant.

No located source promotes the product of 600-cell cover chirality and
chromatic degree to a physical time selector.

## KNOWN

- The all-cover artifact has deterministic SHA-256
  `682e3cfaa0c2912085c0375281817e217f19a54bfc9d6ec9b296844063be7121`.
- It certifies two proper cover orbits, `left` and `right`, exchanged by every
  improper class at the orbit level.
- It certifies that degree sign is proper-even and improper-odd, so `chi=s*d`
  is full-`H4` even.
- The chamber verifier certifies two 60-state `A5` sheets, exchanged by `J`,
  so `gamma` is reflection-odd.
- The full chamber group is `A5 x C2`, while the 600-cell vertex group law is
  `2I`.  Equal cardinality 120 does not identify these carriers.

## Disclosed expected result

The transformation laws already suggest:

1. there are two equivariant correlations `gamma <-> s`, differing by a
   global sign;
2. likewise there are two correlations `gamma <-> d`;
3. there is no equivariant sign-bijection `gamma <-> chi`, because `gamma`
   is reflection-odd and `chi` is reflection-even;
4. symmetry leaves both values of `chi` fixed and therefore selects neither;
5. the regular group laws `2I` and `A5 x C2` cannot be intertwined by a group
   isomorphism.

This is a preregistered negative prediction, not a blind discovery.

## CONTROL

- Pin and parse the all-cover artifact, but independently rebuild the
  icosahedron, its 120 flags, its 60 rotations, central inversion, chamber
  sheets, `gamma` and `J`.
- Enumerate all four maps between every pair of two-point sign sets and test
  equivariance under both rotation and reflection.
- Enumerate the reflection orbits on `(gamma,s)`, `(gamma,d)` and
  `(gamma,chi)` rather than inferring uniqueness from dimensions.
- Reconstruct the `2I` multiplication table from 600-cell quaternions and
  compute its commutator subgroup and involution count.  Compare with the
  exact split chamber group `A5 x C2`.
- Load no Regge action, matter character, Standard-Model chirality, preferred
  cover, phase output, mass or dimensional scale.

## Decision boundary

- `UNIQUE_STATIC_CHIRAL_SELECTOR`: exactly one equivariant construction fixes
  one value of `chi` without extra data.
- `STATIC_CHIRAL_SELECTOR_NO_GO`: the controls pass, `gamma` is odd while
  `chi` is even, both `chi` values remain symmetry-fixed, and any correlations
  with the odd factors come in an unselected pair.
- `OPEN_CONTROL_FAILURE`: any reconstruction or transformation-law control
  fails.

The no-go is scoped to the existing **static** `gamma,J` data.  A new
dynamical term could evade it only if its form and coefficient sign are
derived independently; choosing that sign to favour a desired schedule would
be fitting.  External novelty is **OPEN**.
