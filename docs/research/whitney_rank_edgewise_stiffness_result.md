# Trace stiffness survives the shape-regular control, but no universal factor is selected

Date: 2026-08-11

Preregistration commit: `c7e4335`

Targeted verifier:
`reproducible/verify_whitney_rank_edgewise_stiffness.py`

Targeted result: **9/9 PASS**.  The verifier is registered exactly once.  Only
this targeted calculation was run; the full suite was not run by explicit
user request.

## Headline

The exact Whitney trace-stiffness operator remains well defined and spectrally
controlled on the canonical shape-regular tower

\[
 K_k=\operatorname{Esd}_k(\operatorname{sd}K),
 \qquad k=1,2,4.
\]

The refinement-factor spread across form degrees improves:

\[
 1.26753\longrightarrow1.14892.
\]

Thus the earlier improvement under barycentric refinement was not caused
*only* by degenerating element shapes.

However, the factors do not repeat, and the accumulated degree scales move
farther apart rather than together.

> **PATTERN:** the *increments* become more nearly degree independent on one
> clean shape-regular step.

> **DERIVED NUMERICAL NEGATIVE:** there is no repeated universal refinement
> factor on this control.

No absolute stiffness or physical constant is selected.

## Complete carriers

The control is the closed boundary of the 4-simplex.  Barycentric subdivision
is applied once; edgewise resolution then increases directly.

| edgewise `k` | f-vector | top tetrahedra | duplicated Whitney dimension | ordered element blocks |
|---:|---|---:|---:|---:|
| 1 | `(30,150,240,120)` | 120 | 1,800 | 20 |
| 2 | `(180,1140,1920,960)` | 960 | 14,400 | 53 |
| 4 | `(1320,9000,15360,7680)` | 7,680 | 115,200 | 53 |

The last column counts distinct ordered Gram/integration blocks used by the
implementation, not congruence classes.  Different vertex presentations of a
congruent tetrahedron may produce different ordered blocks.  The geometric
shape theorem remains the three congruence classes certified separately.

All vertices were merged by exact rational barycentric weights.  No geometric
tolerance was used.

## Calibration and solver certificate

The `k=1` carrier is isomorphic to the already certified first barycentric
control.  All six old spectral edges are reproduced with maximum relative
error

\[
 1.44\times10^{-15}.
\]

Across all three carriers and form degrees:

- every shared-face metric agrees exactly;
- every occurrence graph is connected;
- maximum row-image orthonormality residual is `1.11e-15`;
- maximum Ritz residual over all 18 extremal blocks is `8.76e-12`, well below
  the preregistered `1e-7` gate;
- all nine quotient gaps are positive.

The exact constraint row/rank counts are:

| `k` | degree 0 | degree 1 | degree 2 |
|---:|---:|---:|---:|
| 1 | `720 / 450` | `720 / 570` | `240 / 240` |
| 2 | `5760 / 3660` | `5760 / 4620` | `1920 / 1920` |
| 4 | `46080 / 29400` | `46080 / 37080` | `15360 / 15360` |

Each entry is `rows / exact rank`.

## Exact shape scaling and spectral gaps

The stable normalized shape set gives

\[
 \frac{a_4}{a_2}=2
\]

to numerical relative error below `1e-15`, as predicted before the spectral
calculation.  The first jump is larger because `k=2` realizes all three shape
classes whereas `k=1` contains only the parent orthoscheme.

| `k` | worst local Dirac norm `a_k` | gap degree 0 | gap degree 1 | gap degree 2 |
|---:|---:|---:|---:|---:|
| 1 | 11.936870 | 2.293663 | 2.351464 | 2.802834 |
| 2 | 26.233649 | 3.597546 | 3.839674 | 5.572281 |
| 4 | 52.467298 | 6.379754 | 6.732737 | 11.225868 |

At the clean `k=2 -> 4` step, the three gaps grow by approximately

\[
 (1.7734,\;1.7534,\;2.0146).
\]

This is close to the geometric factor two, especially in degree two, but it
is not exact and no exponent is fitted.

## Frozen no-fit comparison

With

\[
 s_{k,p}=a_k/g_{k,p},
\]

the complete scale vectors are

| `k` | degree 0 | degree 1 | degree 2 |
|---:|---:|---:|---:|
| 1 | 5.204283 | 5.076356 | 4.258856 |
| 2 | 7.292096 | 6.832260 | 4.707883 |
| 4 | 8.224031 | 7.792863 | 4.673786 |

Their step ratios are

\[
 R^{1\to2}=(1.401172,\;1.345898,\;1.105434),
\]

\[
 R^{2\to4}=(1.127801,\;1.140598,\;0.992757).
\]

The preregistered factor-repetition diagnostic is

\[
 \frac{R^{2\to4}}{R^{1\to2}}
 =(0.804898,\;0.847462,\;0.898071),
\]

far outside the `1e-6` equality gate.  Exact repetition is refuted.

The ratio-vector spreads are

\[
 \operatorname{spread}(R^{1\to2})=1.267532,
 \qquad
 \operatorname{spread}(R^{2\to4})=1.148919.
\]

That is the precise positive pattern: successive *changes* are more alike
across degrees.

## Hostile reading: the absolute degree mismatch worsens

The scale vectors themselves have degree spreads

\[
 1.22199,\qquad1.54891,\qquad1.75961
\]

at `k=1,2,4`.  Therefore the degrees are not converging to one common value.
More nearly parallel flow does not erase the offsets already accumulated.

This distinction matters:

- the preregistered PATTERN concerns the spread of the **step factors**;
- a stronger claim that one degree-independent stiffness has emerged is
  **false on these data**.

The second-step ratios being close to one are consistent with stabilization
of the dimensionless quantities, but one clean step cannot establish an
asymptotic limit.  That reading remains **OPEN**, not a fitted conclusion.

## Comparison with the degenerating barycentric tower

The old barycentric control gave factor spreads

\[
 3.05498\longrightarrow1.34486.
\]

The new shape-regular control gives

\[
 1.26753\longrightarrow1.14892.
\]

So shape degeneration did strongly contaminate the old magnitudes.  Yet the
direction of improved step-factor balance survives without it.  The honest
conclusion is neither “the old signal was physical” nor “it was entirely a
mesh artefact”: the large old factors were contaminated, while a weaker
finite-level pattern remains.

## Physical status

- **DERIVED:** a complete exact shape-regular control through `k=4`.
- **DERIVED:** exact face agreement, ranks, positive gaps, and factor-two
  local shape scaling.
- **DERIVED NUMERICAL:** all spectral residual gates pass.
- **DERIVED NUMERICAL NEGATIVE:** refinement factors do not repeat.
- **PATTERN:** step-factor spread improves from `1.26753` to `1.14892`.
- **DERIVED NEGATIVE ON THESE LEVELS:** the absolute degree scales do not
  merge; their spread increases.
- **OPEN:** stabilization at further resolutions and on the 600-cell.
- **OPEN:** convergence to a continuum differential operator.
- **OPEN:** absolute stiffness, time, inertia, mass, causal speed, and Planck
  units.

The useful result is methodological and mathematical: the next continuum
question can now be asked without a degenerating mesh.  There is still no new
physical constant.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_rank_edgewise_stiffness.py
```

Expected result: `9/9`.

