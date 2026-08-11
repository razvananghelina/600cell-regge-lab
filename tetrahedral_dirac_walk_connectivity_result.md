# The literal robust walk is unitary but trapped in 1,440 decagons

Date: 2026-08-11  
Preregistration commit: `f52a661`

## Result

The robust three-swap translation transplanted in commit `69447b9` is not a
globally connected dynamics on the barycentric 600-cell chamber carrier.

> **DERIVED CONNECTIVITY NO-GO:** every finite word made from the literal
> robust translation, its inverse and completely arbitrary chamber-local
> coins preserves 1,440 independent blocks.  Each block contains exactly ten
> chambers.

The targeted verifier passes `12/12` in about 0.5 seconds.  No full suite was
run.

## Exact obstruction

The published robust translation crosses chambers only through the two
intrinsic colour involutions

\[
s_2,\qquad s_3.
\]

Their Coxeter product has order five.  The support graph they generate is
therefore a disjoint union of alternating decagons:

\[
14{,}400=1{,}440\times10.
\]

On the full eight-component carrier, each decagon gives an invariant
80-dimensional subspace.  Equivalently, the dynamics has 1,440 exact orbit
projectors.

The corrected macro shift preserves each projector entry by entry.  Its
inverse does too.  An arbitrary chamber-local operator only mixes the eight
components at a fixed chamber, so it also preserves every orbit projector.
Closure under products then proves the result for:

- arbitrary coin angles;
- arbitrary full (8\times8) local mixing;
- position- and time-dependent coins;
- arbitrary mass coins;
- arbitrary finite sequences of forward and inverse translations.

This is a support theorem, not a failed search over coins.

## Geometric control

The full four-coloured (H_4) chamber graph is connected.  Moreover, colours
0 and 1 each cross out of the decagonal blocks at every chamber, and the
resulting quotient graph on the 1,440 blocks is connected.

Therefore the disconnection is not a defect of the 600-cell geometry.  It is
caused precisely by transplanting a flat-lattice translation that uses only
the pair (s_2,s_3).

## Physical meaning

The previous result established a valid local unitary mechanism.  The
present result shows why that was not yet enough:

- **unitarity:** yes;
- **strict finite propagation:** yes;
- **motion through the complete spatial carrier:** no;
- **effective three-dimensional propagation:** no.

The literal published motor runs, but on this geometry it runs around 1,440
separate ten-step tracks.  It cannot be the matter dynamics of the whole
600-cell space.

## Remaining canonical route

The (H_4) Coxeter diagram itself contains exactly three consecutive bonds:

\[
(s_0,s_1),\qquad(s_1,s_2),\qquad(s_2,s_3).
\]

This supplies a target-independent candidate for three ordered spatial
translation stages while using all four chamber directions.  It was
disclosed in the preregistration before the present no-go was computed.

It is not yet a derived physical evolution.  The next gate must test:

1. whether the three substituted robust translations remain exact local
   permutations;
2. whether their periodic support dynamics is connected with the published
   local spin coin;
3. how many order/orientation variants geometry genuinely permits;
4. whether the resulting low-energy/refinement behaviour is Dirac-like,
   without fitting a target dispersion.

## Status ledger

- **DERIVED:** 1,440 causal components, all simple ten-cycles.
- **DERIVED:** 1,440 invariant projectors survive arbitrary local coins.
- **DERIVED NEGATIVE:** the literal two-colour robust transplant is not a
  global matter propagation law.
- **DERIVED CONTROL:** admitting all four intrinsic colours restores chamber
  connectivity.
- **STRUCTURAL:** the ordered three-bond Coxeter schedule is the unique small
  geometry-led continuation currently visible.
- **OPEN:** connectedness of that schedule with a non-fitted coin.
- **OPEN:** relation to Whitney/Kähler--Dirac dynamics and continuum
  refinement.
- **NOT CLAIMED:** a physical particle, a mass derivation, (c), or Planck
  scales.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python \
  reproducible/verify_tetrahedral_dirac_walk_connectivity.py
```

Expected result: `12/12`.
