# Existing chamber chirality cannot select the chromatic `Z2`

Date: 2026-08-17

## Verdict

> **DERIVED EXACT / STRUCTURAL NEGATIVE:** the already derived chamber
> grading `gamma` and reflection `J` do not select either value of the
> residual chromatic invariant
>
> ```text
> chi = cover chirality * chromatic-degree sign.
> ```

The mechanical outcome is

```text
STATIC_CHIRAL_SELECTOR_NO_GO
```

The targeted verifier passed `14/14` twice.  The full repository suite was
not run.

## Provenance

- prior-art gate: `e946b0a`;
- preregistration: `70525f1`;
- verifier registered before first execution: `17deed8`;
- verifier:
  `reproducible/verify_gravity_600cell_chromatic_chamber_selector.py`;
- result:
  `reproducible/gravity_600cell_chromatic_chamber_selector.json`;
- deterministic result SHA-256:
  `ad0094cec72ac10c7f5d5e098b1d5f11ce3743cc83e3d73013744d04c6a4ea5c`.

The negative transformation-law prediction was disclosed before the test.
No Regge output, matter target or preferred chirality was loaded.

## The decisive mismatch

The four signs transform under reflection as

```text
gamma -> -gamma,
s     -> -s,       s = left/right cover chirality,
d     -> -d,       d = chromatic-degree sign,
chi=s*d -> chi.
```

Therefore `gamma`, `s` and `d` are reflection-odd, while `chi` is
reflection-even.  Exhausting all four functions between two sign sets gives:

```text
equivariant bijections gamma -> s    : 2
equivariant bijections gamma -> d    : 2
equivariant bijections gamma -> chi  : 0
```

The first two counts are nonselection: one may correlate the signs equally
well or oppositely.  The two choices differ by the global sign that the
geometry does not fix.  The last count is a genuine parity obstruction: an
odd sign cannot be bijectively identified with an even sign.

Direct orbit enumeration says the same thing.  Each of

```text
(gamma,s), (gamma,d), (gamma,chi)
```

has two reflection orbits.  Moreover `chi=+1` and `chi=-1` are two separate
singleton orbits because reflection fixes each.  An invariant potential may
assign different values to those two orbits, but symmetry supplies neither
the values nor which one is lower.  Choosing that coefficient sign after
inspecting dynamics would be fitting.

## The two 120-state carriers are not the same group object

The verifier rebuilt the chamber geometry without parsing its old result:

```text
icosahedron f-vector       (12,30,20)
complete flags             120
proper rotations            60
orientation sheets       60+60
full chamber action         120, free and transitive
```

It then reconstructed both order-120 multiplication laws exactly enough to
compute discrete invariants:

| regular group law | nonidentity involutions | commutator subgroup |
|---|---:|---:|
| binary icosahedral `2I` | 1 | 120 |
| chamber `A5 x C2` | 31 | 60 |

Thus `2I` is perfect and non-split, while the chamber group has a nontrivial
determinant quotient.  They are not isomorphic.  A dimension equality

```text
120 vertices = 120 chambers
```

cannot supply a group-law-preserving identification.  Arbitrary vector-space
bijections exist, but using one to correlate chirality would add precisely
the basepoint/basis choice excluded by the protocol.

This statement is deliberately scoped: extra geometric data could define a
different correspondence.  None is currently derived.

## Hostile framing audit

### Could `J` itself choose the sign?

No.  `J` is the involution exchanging the two chamber sheets; it does not
prefer either sheet.  On the chromatic side the corresponding improper
symmetry swaps `s` and `d` simultaneously and leaves `chi` unchanged.

### Could one use `gamma*s` as a reflection-even coupling?

Yes as an algebraic term, but its two invariant values are both allowed.
Neither symmetry nor the existing static geometry fixes the sign of its
coefficient.  The same applies to `gamma*d`.  This is an admissible future
dynamical question, not a current derivation.

### Does observed weak chirality choose it?

Not without a separately derived map from the matter representation to this
cover invariant.  Selecting the sign from the Standard Model would be target
fitting, especially because the orbifold incidence matter route has already
failed its canonical-kernel test.

## Prior-art status

The non-split `2.A5 ~= SL(2,5)` extension, the split full icosahedral group and
the homogeneous-set map criterion are standard group theory.  Fisk's ten
colourings are in [Coloring the 600
Cell](https://arxiv.org/abs/0802.2533).  The closest internal theorem,
`hopf_axis_orientation_verdict.md`, had already proved an analogous static
orientation no-go for handed Hopf fibres.

The post-result search found no primary source assigning physical time
chirality to the specific invariant `chi`.  This computation is chiefly a
repository consistency/no-fitting result; external novelty remains **OPEN**
and is not claimed.

## Status ledger

- **DERIVED EXACT:** `gamma`, cover side `s` and degree sign `d` are
  reflection-odd; `chi=s*d` is even.
- **DERIVED EXACT:** two equivariant `gamma<->s` correlations and two
  `gamma<->d` correlations remain.
- **DERIVED EXACT:** zero equivariant sign-bijections `gamma<->chi`.
- **DERIVED EXACT:** both `chi` values are separately symmetry-fixed.
- **DERIVED EXACT:** the regular group laws `2I` and `A5 x C2` are
  nonisomorphic.
- **STRUCTURAL NEGATIVE:** existing static chamber data do not select the
  chromatic chiral class.
- **OPEN:** a coefficient-free dynamical parity-breaking term.
- **OPEN:** any canonical additional bridge between the two carriers.
- **OPEN:** external novelty.

## Programme consequence

The cheap escape route is closed.  We cannot repair the schedule ambiguity by
declaring that the old chamber `gamma` chooses it.

The next technically defensible route is triangulation independence.  The
observed quadratic schedule defect is exactly the kind of broken discrete
gauge symmetry for which improved/perfect Regge actions were developed.  The
primary references are Bahr and Dittrich, [Improved and Perfect Actions in
Discrete Gravity](https://arxiv.org/abs/0907.4323), and Dittrich and
Steinhaus, [Path integral measure and triangulation independence in discrete
gravity](https://arxiv.org/abs/1110.6866).

That route is harder and may require nonlocal boundary terms in four
dimensions.  But it attacks the actual failure—schedule dependence—instead
of choosing one side of it.
