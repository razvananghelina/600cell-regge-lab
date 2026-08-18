# Chromatic compatibility leaves an exact chiral `Z2`

Date: 2026-08-17

## Verdict

> **DERIVED COMPUTATIONAL / STRUCTURAL NEGATIVE:** all ten five-colourings of
> the 600-cell have chromatic degree of magnitude `72`, but the
> orientation-compatible ordered colourings do **not** form one canonical
> class.  They split under the 7200 proper rotations into two orbits of size
> `300`, supported respectively on the five left-coset and five right-coset
> covers.

The frozen mechanical outcome is

```text
CHIRAL_COVER_AMBIGUITY
```

Thus the chromatic degree removes the arbitrary even/odd phase ordering for a
*fixed* cover, but it moves the remaining ambiguity to a left/right chiral
choice of cover.  It does not yet select a unique nonlinear tick.

## Provenance and reproducibility

- prior-art gate: `1071c40`;
- preregistration: `a1d9bf0`;
- verifier registered before first execution: `e6b3a0a`;
- removal of volatile runtime fields only: `3118380`;
- post-result explicit `Z2` diagnostic, without changing the outcome rule:
  `5aa5539`;
- verifier:
  `reproducible/verify_gravity_600cell_chromatic_cover_orbits.py`;
- result:
  `reproducible/gravity_600cell_chromatic_cover_orbits.json`;
- deterministic result SHA-256:
  `682e3cfaa0c2912085c0375281817e217f19a54bfc9d6ec9b296844063be7121`.

The final targeted verifier passed `22/22` twice, with the same artifact hash.
The full repository suite was not run.

The prediction of two 300-element orbits was disclosed before the
calculation.  The result is therefore a confirmed prediction, not a blind
discovery.

## Exhaustive input rather than ten convenient examples

The verifier did not assume Fisk's list of ten colourings.  It rebuilt the
full order-14400 `H4` action, then obtained:

```text
25 maximum independent sets of size 24;
alpha(G) = 24                 [CP-SAT status OPTIMAL];
no further size-24 set        [CP-SAT status INFEASIBLE];
10 exact covers by five cells;
5 left-coset + 5 right-coset covers.
```

Because every class in a proper five-colouring has size at most 24 and the
five classes contain all 120 vertices, every class must have size exactly 24.
The exclusion calculation therefore makes the ten-cover census exhaustive.

Independently, every generated symmetry was reconstructed as a real `4 x 4`
orthogonal matrix.  Its determinant produced exactly

```text
7200 proper rotations;
7200 improper symmetries.
```

The orientation classification agreed with, but was not defined by, the
presence or absence of quaternion conjugation.

## All ten exact degrees

For every cover and all `5!` orders, the verifier pushed the complete signed
600-tetrahedron fundamental chain to each of the five target facets.  All five
degree candidates agreed with all five direct signed-preimage counts.

The complete multiset is

```text
degree -72 : 600 ordered colourings;
degree +72 : 600 ordered colourings;
degree   0 :   0 ordered colourings.
```

Every unordered cover has degree magnitude `72`; its 120 orders divide into
60 of each sign.  The number `72` is therefore not an accident of the first
left-coset cover.

## Exact orbit census

On unordered covers:

```text
proper H4+ orbits : 5 left, 5 right;
full H4 orbit     : all 10;
improper map      : left <-> right.
```

On positive-degree ordered covers:

```text
proper H4+ orbits : 300 + 300;
stabilizer        : 24 for each;
time reversal     : preserves each orbit.
```

The same `300+300` split holds in the negative-degree sector.  On all 1200
orders, full `H4` still has two 600-element orbits.

There is also a direct algebraic explanation.  A proper quaternionic symmetry
has the form

```text
R_(l,r)(q) = l q r^-1.
```

It maps a left-coset cover to a left-coset cover because

```text
l (g H) r^-1 = (l g r^-1) (r H r^-1),
```

and analogously maps a right-coset cover to a right-coset cover.  Quaternion
conjugation exchanges left and right.  At the same time, proper maps preserve
the source orientation and chromatic degree, whereas conjugation has
determinant `-1` and reverses the degree.

If left/right cover chirality is encoded by `s=+1/-1` and the degree sign by
`d=+1/-1`, then

```text
chi = s*d
```

is invariant under the **full** `H4`: a reflection flips both factors.  The
two 600-element full-group orbits are exactly `chi=+1` and `chi=-1`.  This is
why simply admitting reflections does not collapse the two classes.

## Hostile framing audit

### Could one quotient by reflections and declare one class?

No, not while retaining the positive-degree rule on a fixed oriented carrier.
Improper maps send positive degree to negative degree, so the positive sector
is not closed under that quotient.  If the spatial orientation is transported
as additional data, the invariant correlation `chi` still leaves two full-
group orbits.  Identifying them would require an extra odd relabelling of the
phase order.  That relabelling changes the staircase triangulation and is not
currently a derived symmetry.

### Is one chirality physically forbidden?

**OPEN.**  Pure Regge gravity with dust supplies no parity-violating term that
chooses left over right.  The repository has chiral structures (`gamma`, `J`
and matter representations), but no derived map currently identifies their
chirality with `chi`.  Using the Standard Model's observed chirality to pick a
cover sign before deriving such a map would be fitting.

### Does this refute the integer 72?

No.  The magnitude `72` is exact and universal over all ten covers.  What is
refuted is the stronger claim that this invariant alone selects one global
schedule class.

## Post-result prior-art audit

Fisk proves the ten five-colourings in [Coloring the 600
Cell](https://arxiv.org/abs/0802.2533).  The left/right isoclinic distinction
of embedded 24-cells is also described in modern treatments of the 600-cell,
including the explicit coset construction in [Geometry and combinatorics of
the 600-cell and the
120-cell](https://repository.tudelft.nl/file/File_fd61d63c-22a4-46f7-8992-cf7b511df139).
The degree formalism itself is standard in [Toric residue and combinatorial
degree](https://arxiv.org/abs/math/0309409).

The post-result search did not locate the explicit universal magnitude `72`,
the `300+300` positive-degree orbit census, or the invariant `chi=s*d` in a
primary source.  Search cannot prove novelty; their external novelty remains
**OPEN**.

## Status ledger

- **DERIVED COMPUTATIONAL:** exhaustive `25`-cell, ten-cover census.
- **DERIVED EXACT:** all ten chromatic degree magnitudes equal `72`.
- **DERIVED EXACT:** the complete degree multiset is `600*(-72)+600*(+72)`.
- **DERIVED COMPUTATIONAL:** proper rotations split covers as five left plus
  five right and compatible schedules as `300+300`.
- **DERIVED / STRUCTURAL:** `chi = cover chirality * degree sign` is the exact
  full-`H4` `Z2` invariant.
- **STRUCTURAL NEGATIVE:** chromatic orientation alone does not select one
  canonical ordered cover.
- **OPEN:** whether an independently derived chiral physical structure fixes
  `chi`.
- **OPEN:** physical force of chromatic compatibility itself.
- **OPEN:** external novelty of `72`, `300+300` and the explicit `Z2`.

## Consequence for the dynamics programme

The bare order-24 Regge+dust quotient already breaks nonlinear schedule
covariance at quadratic order.  Product orientation does not choose a
schedule, and chromatic orientation now leaves a chiral `Z2`.  Therefore the
current geometry still does not define a unique nonlinear tick.

The next defensible branch is not to choose the favourable chirality.  It is
either:

1. derive a target-blind coupling from an already existing chiral operator
   (`gamma` or `J`) to `chi`; or
2. seek a schedule-independent/perfect Regge action whose boundary evolution
   is invariant under the local staircase flips.

Without one of these, claims about a derived `c`, Planck time, Planck mass or
particle masses remain premature.

### Subsequent closure of branch 1

The targeted follow-up
`gravity_600cell_chromatic_chamber_selector_result.md` closes the first cheap
option for the *existing static* chamber data.  `gamma` is reflection-odd,
whereas `chi` is reflection-even; there are zero equivariant sign-bijections
between them, and both `chi` values remain fixed.  Consequently `gamma,J`
alone do not select the residual class.  A genuinely new derived dynamical
coupling would be required.  The schedule-independent/perfect-action branch
is therefore the next active route.
