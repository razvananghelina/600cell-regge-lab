# The five-colouring has exact chromatic degree 72

Date: 2026-08-17

## Verdict

> **DERIVED MATHEMATICAL / STRUCTURAL:** the fixed five-colouring of the
> oriented 600-cell defines a nonzero map
> `K -> boundary(Delta^4)` of degree `-72` in the preregistered label
> convention.  Relabelling by an odd permutation reverses the degree.  Thus
> the two 60-element staircase orbits are exactly the two orientations of the
> chromatic four-simplex.

Frozen outcome:

```text
CHROMATIC_ORIENTATION_LINE_DERIVED
```

The targeted verifier passed `10/10`.  The full suite was not run.

### Subsequent all-cover scope correction

The exhaustive follow-up
`gravity_600cell_chromatic_cover_orbits_result.md` computed all ten
five-colourings.  Every one has degree magnitude `72`, but the compatible
ordered covers split under proper rotations into two chiral orbits of size
`300`, carried by the five left- and five right-coset covers.  Thus the fixed-
cover orientation line is exact, while its stronger use as a globally unique
schedule selector is **REFUTED / STRUCTURAL NEGATIVE**.

This is not yet a physical schedule selector.  No current Regge, causal or
matter axiom requires the staircase order to agree with the induced chromatic
orientation.

## Provenance

- prior-art gate: `38cce14`
- preregistered protocol: `4458c23`
- verifier registered before first execution: `92dd8f8`
- verifier:
  `reproducible/verify_gravity_600cell_chromatic_degree_selector.py`
- result:
  `reproducible/gravity_600cell_chromatic_degree_selector.json`
- result SHA-256:
  `a3847065571cecc37fb7f5a68896287c5c900c572d61d6871eac8594b5556332`

No Regge action, nonlinear output, continuum target or preferred sign was
loaded.

## Exact degree calculation

The determinant-oriented source chain is closed exactly.  Its minimum
absolute tetrahedral determinant is

```text
0.15450849719465062,
```

against the preregistered ambiguity threshold `1e-10`.

For target facet `i`, the pushforward coefficients are

```text
P = (-72, +72, -72, +72, -72).
```

The oriented boundary of the target four-simplex has coefficients

```text
(+1, -1, +1, -1, +1).
```

Therefore all five independent degree candidates are

```text
((-1)^i P_i) = (-72, -72, -72, -72, -72).
```

Each target facet has exactly 120 source tetrahedral preimages.  On every
facet, 96 local degrees have sign `-1` and 24 have sign `+1`, giving the same
signed count

```text
24 - 96 = -72.
```

This independent preimage count matches the chain pushforward on all five
facets.

## Complete relabelling and symmetry census

All 120 total colour orders obey exactly

```text
degree(order) = sign(order) * degree(identity).
```

Consequently:

```text
60 even orders : degree -72
60 odd orders  : degree +72.
```

The rebuilt setwise `H4` cover action again contains 1440 actions, induces
the 60 even permutations and has kernel 24.  All 60 induced `A5`
permutations preserve degree `-72` exactly.

Thus the earlier two schedule orbits have a precise topological meaning: they
are not the two orientations of the spacetime product chain, but they are the
two orientations of the target chromatic four-simplex.

## What is and is not selected

Given the oriented 600-cell and the *unoriented* five-colouring, the nonzero
degree canonically induces an orientation on the abstract colour simplex: call
the target orientation positive when the colouring map has positive degree.
This reduces the 120 total orders to one 60-element `A5` orbit of
orientation-compatible orders.  No experimental target or Regge outcome is
needed for that mathematical convention.

However, a staircase triangulation of `K x I` is valid for either orbit and
has the same relative fundamental chain.  Standard prism topology supplies no
condition saying that the temporal phase order must orient the auxiliary
colour simplex in the induced way.  Imposing that compatibility is an extra
canonicity axiom unless it can be derived from causality, the action or matter.

Therefore:

- **DERIVED:** a target-blind chromatic orientation line exists;
- **STRUCTURAL:** it distinguishes exactly the two nonlinear schedule
  parities;
- **OPEN:** whether physics selects the compatible rather than incompatible
  chromatic orientation;
- **FORBIDDEN FITTING:** declaring the sign physical after inspecting which
  schedule gives a preferred dynamical answer.

## Prior-art audit

The mathematical ingredients are known.  Fisk derives the ten five-colourings
of the 600-cell in [Coloring the 600
Cell](https://arxiv.org/abs/0802.2533).  Soprounov defines the combinatorial
degree of a colouring map to a simplex boundary and proves its alternation in
the colour order in [Toric residue and combinatorial
degree](https://arxiv.org/abs/math/0309409).  The post-result search did not
locate the explicit value `72` for this 600-cell colouring.  Search cannot
prove novelty, so external novelty of that integer remains **OPEN**.

## Status ledger

- **DERIVED EXACT:** chromatic degree magnitude `72`.
- **DERIVED EXACT:** five facet calculations agree and each decomposes as 24
  positive versus 96 negative local preimages.
- **DERIVED EXACT:** the 120 orders split as `60 x (-72)` and `60 x (+72)`.
- **DERIVED EXACT:** the fixed-cover `A5` action preserves one degree sign.
- **STRUCTURAL:** the chromatic orientation, unlike the product orientation,
  distinguishes the two schedule orbits.
- **OPEN:** physical force of that distinction.
- **DERIVED EXACT (subsequent audit):** all ten colourings have degree
  magnitude `72`.
- **STRUCTURAL NEGATIVE (subsequent audit):** chromatic-compatible schedules
  form two proper-rotation orbits of size `300`; the ambiguity moves to the
  left/right chirality of the cover.

## Next decisive calculation

The all-cover census is complete.  The next admissible question is whether an
independently derived chiral operator fixes the residual `Z2`, or whether a
schedule-independent action removes the need to choose it.
