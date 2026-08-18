# Three-dimensional Whitney flux is not selected by the current axioms

Date: 2026-08-11

Preregistration commit: `64af02c`

Targeted verifier:
`reproducible/verify_whitney_3d_flux_canonicity.py`

Targeted result: **8/8 PASS**.  The verifier is registered exactly once.  The
full suite was not run, by explicit user request.

## Headline

The successful circle flux does not extend uniquely to the
three-dimensional shape-regular carrier using the presently supplied
geometric conditions.

> **DERIVED NEGATIVE FOR UNIQUENESS:** equal counting recovery and diagonal
> Whitney-metric recovery are distinct exact, positive,
> strict-occurrence-local, natural left inverses on `Esd_2(sd K)`.

> **DERIVED LOCALITY CONFLICT:** exact Hilbert adjointness does select a unique
> recovery, but that recovery has an exactly nonzero coefficient outside the
> occurrence star and numerical support reaching the diameter of the refined
> control complex.

Therefore no three-dimensional flux spectrum was computed.  Choosing the
candidate with the more attractive spectrum would be fitting after the
ambiguity is known.

## Complete setting

The controls are the closed boundary of a 4-simplex, subdivided once
barycentrically to produce its canonical face-rank order, then directly
edgewise refined at resolutions `k=1,2`.  They use the exact local Euclidean
Whitney masses of the already certified shape-regular tower.

For degree `p`, `J_p` copies every global cochain coefficient into all
incident tetrahedra.  A strict-occurrence-local recovery has the form

\[
 (L_px)_s=\sum_{T\supset s}w_{T,s}x_{T,s},
 \qquad
 \sum_{T\supset s}w_{T,s}=1.
\]

Orientation signs are included explicitly in the matrices.  The exact sums,
not floating residuals, certify `L_pJ_p=I`.

The two local rules frozen before evaluation were

\[
 w^C_{T,s}={1\over |O(s)|}
\]

and

\[
 w^D_{T,s}={q_{T,s}\over\sum_{U\supset s}q_{U,s}},
 \qquad
 q_{T,s}=(M^{T}_{p})_{s,s}>0.
\]

The second rule is precisely the Hilbert-adjoint recovery for the diagonal,
mass-lumped Whitney metric.  Thus it is not an arbitrary list of weights:
it is constructed functorially from the distinguished oriented-simplex basis
and the supplied local metric.

## Exact local ambiguity

The carrier controls are

| `k` | f-vector | tetrahedra | exact element-metric types |
|---:|---:|---:|---:|
| 1 | `(30,150,240,120)` | 120 | 20 |
| 2 | `(180,1140,1920,960)` | 960 | 53 |

Both recoveries are exact left inverses in every degree and every local mass
diagonal is strictly positive.

| `k` | degree | global/local dimensions | occurrence range | rows with `L^C != L^D` | differing coefficients | maximum exact absolute difference |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0 | 30 / 480 | 12--24 | 0 | 0 | 0 |
| 1 | 1 | 150 / 720 | 4--6 | 0 | 0 | 0 |
| 1 | 2 | 240 / 480 | 2 | 0 | 0 | 0 |
| 1 | 3 | 120 / 120 | 1 | 0 | 0 | 0 |
| 2 | 0 | 180 / 3840 | 12--24 | 0 | 0 | 0 |
| 2 | 1 | 1140 / 5760 | 4--6 | **480** | **2640** | **1/15** |
| 2 | 2 | 1920 / 3840 | 2 | **960** | **1920** | **5/32** |
| 2 | 3 | 960 / 960 | 1 | 0 | 0 | 0 |

The first exact edge witness has six copies.  Counting gives `1/6`; one
local diagonal mass is `1/60`, for which diagonal recovery gives `1/9`.
Their difference is exactly `-1/18`.

The first exact face witness has two copies.  Counting gives `1/2`; a local
diagonal mass `6/5` gives diagonal weight `9/22`, with difference `-1/11`.

Agreement at `k=1` was therefore only a symmetry accident of that level, not
a theorem.  The first nontrivial edgewise refinement exposes the distinction.
Degree 3 remains the required one-copy identity control.

## Exact metric adjoint: unique but nonlocal

For the complete local Whitney metric, adjointness uniquely fixes

\[
 L^A_p=(J_p^*M_{\rm loc,p}J_p)^{-1}J_p^*M_{\rm loc,p}.
\]

All solves and all numerical left-inverse identities have maximum residual
below `2.64e-15`, far inside the preregistered `1e-10` gate.  With relative
support threshold `1e-11`, the observed support is

| `k` | degree | significant coefficients off occurrence support | reached graph radius / graph diameter |
|---:|---:|---:|---:|
| 1 | 0 | 13,440 | 2 / 3 |
| 1 | 1 | 103,680 | 4 / 5 |
| 1 | 2 | 103,680 | 9 / 10 |
| 1 | 3 | 0 | 0 / 0 |
| 2 | 0 | 686,880 | 6 / 6 |
| 2 | 1 | 6,544,800 | 9 / 9 |
| 2 | 2 | 6,949,440 | 19 / 19 |
| 2 | 3 | 0 | 0 / 0 |

The support conclusion is not merely a tolerance artefact.  After the frozen
numerical audit identified a robust coefficient, it was reconstructed over
the exact rational matrices.  At `k=1`, degree zero, row 5 of the adjoint has
coefficient

\[
 {243\over7480}\ne0
\]

on local copy 52, whose underlying global vertex is 3 rather than 5.  The
floating value agrees with the exact rational to `1.39e-17`.

This exact witness was a post-result strengthening of the preregistered
numerical locality test.  Selecting a witness after seeing a matrix is valid
for a falsifying existence proof; it is not used as statistical evidence or
as a fitted physical parameter.

## Attack on the word “canonical”

There are two distinct meanings:

1. a formula is natural and contains no arbitrary labels or fitted numbers;
2. the axioms uniquely select that formula over every other natural formula.

Both `L^C` and `L^D` satisfy the first meaning.  Their exact disagreement
refutes the second under the frozen hypotheses.  In fact, once their local
norms differ, normalized positive integer powers of those norms provide still
more coefficient-free natural rules.  Counting and diagonal weighting are
enough to establish the logical failure; they are not claimed to enumerate
all possibilities.

Adding **exact full-metric adjointness** resolves uniqueness, but selects
`L^A`, whose inverse assembled mass makes it nonlocal.  Hence the current
inputs give a sharp fork:

- strict local tick plus naturality: more than one recovery;
- exact Whitney adjointness: one recovery, but global support.

This does not prove that no additional physical axiom could select a local
map.  Possible new requirements include a commuting-cochain property,
conservation law, or a variational principle formulated without the global
inverse.  None is currently derived by the theory, and one may not be chosen
because of its eventual spectrum.

## Consequence for the circle positive

The circle result remains correct.  There every vertex has exactly two
geometrically exchanged copies; exchange symmetry and reproduction force
weights `(1/2,1/2)`.  The three-dimensional refined carrier has inequivalent
incident tetrahedral extensions, so the same proof no longer applies.

Thus the circle established that finite stiffness *can* be repaired by a
consistency flux.  It did not establish that the 3D geometry selects which
flux.

## Status ledger

- **DERIVED:** both frozen local rules are positive natural exact left
  inverses.
- **DERIVED NEGATIVE:** they differ in degrees 1 and 2 at `k=2`.
- **DERIVED:** the exact full-metric adjoint is the unique metric-adjoint left
  inverse.
- **DERIVED:** one exact rational coefficient proves it is not
  strict-occurrence-local.
- **DERIVED NUMERICAL:** its support reaches the complete simplex-graph
  diameter at `k=2` in degrees 0, 1 and 2.
- **STRUCTURAL NEGATIVE:** the present axioms do not select a bounded-star 3D
  flux.
- **OPEN:** whether a new independently motivated cochain or conservation
  axiom uniquely selects one.
- **NOT CLAIMED:** failure of every conceivable 3D discretization or of
  Kähler--Dirac theory itself.
- **NOT CLAIMED:** physical time, inertia, mass, causal speed or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_3d_flux_canonicity.py
```

Expected result: `8/8`.

