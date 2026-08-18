# SUPERSEDED STEP 3: comparison of an incomplete operator list

> **DERIVED CORRECTION (2026-08-09):** this file's `N=3` premise is false.
> The twisted-kernel formula averaged over whole stabilizers and artificially
> cancelled every nontrivial character.  The corrected census has 62 maps and
> one short complex.  See `orbifold_incidence_route.md` and
> `reproducible/verify_orbifold_incidence_route.py`.  This historical text is
> retained so the failed preregistration remains auditable.

Preregistration commit: `36bd6825b6b64edc6c95c1ac67a21a80693318fc`
That commit contains the full operator list and its kernel/cokernel characters,
and states in its message that no target comparison had been performed.  The
enumeration script never reads a target, so it could not have branched on one.

## Result: no canonical incidence operator carries the target

`N = 3` canonical nonzero incidence operators exist.  Their kernels and
cokernels, in irrep order with dimensions `[1,2,2,3,3,4,4,5,6]`:

| operator | rank | ker | coker |
|---|---:|---|---|
| edge -> vertex | 12 | `2x4' + 2x5` (dim 18) | `0` |
| face -> vertex | 12 | `2x4'` (dim 8) | `0` |
| face -> edge | 20 | `0` | `2x5` (dim 10) |

The target is `M16 = 2 rho_1 + rho_2 + 2 rho_5 + rho_9`, i.e. multiplicities on
the irreps of dimension `1, 2, 3, 6`.  All three same-dimension labelling
variants were tested, together with `3 x M16`.

**Matches: 0 out of 3 operators x 3 labelling variants x 2 (kernel, cokernel).**

The failure is not marginal.  The kernels are supported entirely on the
four- and five-dimensional irreps, while `M16` is supported entirely on the
one-, two-, three- and six-dimensional ones.  The supports are disjoint.

## The three-term complex does not exist either

For a middle cohomology one needs `d1 . d2 = 0`.  Exactly:

`d1 . d2 != 0`, with `max |entry| = 1920`.

Rescaling cannot repair this: the composition is bilinear in the two scalars,
so it vanishes only if it already does.  Hence there is no canonical
`face -> edge -> vertex` complex and no middle cohomology to compute.
For the record, `dim ker d1 = 18`, `dim im d2 = 20`, `dim(ker ∩ im) = 8`.

## Verdict

**KILL BOUNDARY REACHED for the canonical construction.**  The complete family
of incidence operators between induced line modules from the derived
stabilizers `C10`, `C4`, `C6`, supported on the incidence double coset, cannot
produce `M16`, three copies of it, or any complex whose cohomology could.

**DERIVED NEGATIVE**, with this scope: line modules only; the three derived
stabilizers only; operators supported on the incidence double coset only.  Not
covered, and therefore still open in principle: higher-rank equivariant
bundles, other subgroups, and any different notion of canonicity.  Operators
supported off the incidence coset are excluded by definition, not by evidence.

Note that the twisted operators did not have to be argued away: for every
non-trivial character pair the incidence sum vanishes identically, which is why
`N = 3` rather than the naive `10*4 + 10*6 + 4*6`.

## What this closes

`math_to_physics_bridge_audit.md` called the full `(2,3,5)` orbifold "the
surviving route" after the `C10`-only functor was refuted.  With the canonical
reading of "incidence operator" it is now closed as well.  The repository has a
discrete geometric precursor and no derived matter functor.
