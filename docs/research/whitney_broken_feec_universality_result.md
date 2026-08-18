# The projection ambiguity survives refinement and worsens in top degree

Date: 2026-08-11

Preregistration commits: `72d4eaa`, with pre-result gate clarification
`230195b`.

Targeted verifier:
`reproducible/verify_whitney_broken_feec_universality.py`

Targeted result: **8/8 PASS**.  The verifier is registered exactly once.  The
full suite was not run, by explicit user request.

## Headline

Removing the weak-penalty spurious modes does not make the two
three-dimensional broken-FEEC projections converge toward one another on the
available refinement step.

> **DERIVED NEGATIVE FOR UNIFORM MICROSCOPIC CONVERGENCE:** the exact maximum
> local weight differences remain `1/15` and `5/32` from `k=2` to `k=4`, while
> the fractions of affected simplices increase.

> **PATTERN NEGATIVE FOR A COMMON LOW-ENERGY FLOW:** in the coefficient-free
> infinite-penalty constrained spectrum, the top-form low quartet separates
> by `0.1048%` at `k=2` but by `1.4378%` at `k=4`, a factor `13.72` larger.

This does not prove distinct continuum limits from two levels.  It does remove
the only current numerical support for saying that the projection ambiguity
is becoming irrelevant.

## Why this spectrum is the right diagnostic

The preceding `alpha=1` control had the correct kernel but lay outside the
published spectral-convergence regime.  The broken-FEEC theorem assumes

\[
 \alpha_h\geq C h^{-s},
\]

and the authors explicitly observe spurious low weak-penalty eigenvalues.

This audit instead takes the exact strong-penalty limit at each fixed mesh.
Finite eigenbranches then lie in the conforming image `im(J)`, and the pencil
contains no stabilization coefficient:

\[
 \widehat K^X_{h,p}=J_p^T K^X_{h,0,p}J_p,
 \qquad
 \widehat M_{h,p}=J_p^TM_{h,p}J_p,
 \qquad X\in\{C,D\}.
\]

Thus no arbitrary `alpha`, branch threshold, target eigenvalue or fitted scale
enters the comparison.

This is still an effective constrained limit, not a uniquely selected finite
microscopic tick.

## Exact microscopic result

The two local recoveries are equal counting (`C`) and diagonal Whitney
mass-lumped (`D`).  Their exact comparison is

| degree | `k` | global rows | differing rows | row fraction | differing coefficients | exact maximum weight difference |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 2 | 1,140 | 480 | 42.105% | 2,640 | `1/15` |
| 1 | 4 | 9,000 | 4,800 | 53.333% | 26,400 | `1/15` |
| 2 | 2 | 1,920 | 960 | 50.000% | 1,920 | `5/32` |
| 2 | 4 | 15,360 | 11,520 | 75.000% | 23,040 | `5/32` |

The first exact witnesses also repeat unchanged:

- degree 1: counting `1/6`, diagonal `1/9`, difference `-1/18`;
- degree 2: counting `1/2`, diagonal `9/22`, difference `-1/11`.

Because the exact suprema do not decrease, the projection matrices do not
approach one another coefficientwise in supremum norm on the complete broken
carrier.  The growing row fractions show that the obstruction is not confined
to a fixed exceptional set.

This statement is exact on the two frozen levels.  An all-level repetition
theorem is not claimed.

## Strong-limit low spectrum

All controls retain Betti nullities

\[
 (1,0,0,1)
\]

for both candidates.  The maximum relative Ritz residual is `8.88e-13`, more
than five orders of magnitude inside the preregistered `1e-7` gate.

### Degrees zero and one

The constrained pencils agree exactly on the controls.  In degree zero there
is no lower coderivative term.  In degree one that term depends only on the
degree-zero projection, for which the two recovery rules agree.  These are
positive controls, not evidence that every degree is universal.

### Degree two

The constrained matrices are different—at `k=4` their difference has
4,396,440 significant entries—but the first four eigenvalues agree to maximum
relative discrepancy `3.72e-14`:

\[
 \lambda^{C}_{1,\ldots,4}=2.914823509127\ldots,
 \qquad
 \lambda^{D}_{1,\ldots,4}=2.914823509127\ldots.
\]

Higher reported modes already differ (`7.2950...` versus `7.6454...`).  The
low-quartet agreement is therefore a symmetry-protected or accidental
isospectral sector on this control, not equality of operators.

### Degree three: the preregistered failure

The first positive eigenvalue is fourfold degenerate in both candidates.

| `k` | counting | diagonal Whitney | relative difference |
|---:|---:|---:|---:|
| 2 | 4.374103306795 | 4.369517911725 | 0.10483% |
| 4 | 4.334310151821 | 4.271991258800 | 1.43780% |

Every member of the first quartet has the same relative separation to the
shown precision.  The ratio of the two relative differences is

\[
 13.7155.
\]

The preregistered common-flow gate therefore fails in degree three.

## Interpretation

Three logically separate statements now coexist:

1. **Topology:** both projections preserve the exact `S^3` harmonic kernel.
2. **Microscopic operator:** their local coefficients remain separated by
   fixed rational amounts under the tested refinement.
3. **Low-energy constrained operator:** some sectors agree, but the top-form
   sector separates more strongly at the new level.

The first statement is a real robustness result.  It does not imply the other
two.

The present evidence therefore does **not** justify saying “the choice of
projection washes out in the continuum.”  To establish that claim would now
require an analytic approximation theorem for both `P` and `P*` on this
specific simplicial manifold, capable of explaining why the observed
top-degree separation must eventually reverse.  No such proof currently
exists in the repository.

## Attack on the negative framing

One worsening step is not a no-go theorem for a common continuum limit.
Finite-size effects, the piecewise-flat singular geometry, or a later reversal
could still occur.  Fitting a divergence exponent from `k=2,4` would be just as
invalid as fitting convergence.

Nor does coefficientwise nonconvergence on the complete broken space exclude
low-energy universality after mismatch modes are removed.  That is why the
exact microscopic negative and the spectral PATTERN negative are reported
separately.

Finally, the limiting geometry here is the fixed piecewise-flat boundary of a
4-simplex after subdivision, not the round three-sphere.  No round-`S^3`
eigenvalue comparison is being smuggled into the test.

## Physical consequence

The duplicated Whitney construction now has a sound, known mechanism for
preserving topology, but still lacks a derived dynamics:

- the finite conforming projection is not unique;
- the ambiguity does not shrink microscopically on the tested tower;
- its top-degree effect grows rather than shrinks in the strong low-energy
  control;
- the strong constraint itself is an effective limit, not a microscopic
  causal update.

Thus there is still no route here to a uniquely derived speed, inertia, mass,
time quantum or Planck scale.  The correct next mathematical question is no
longer another unconstrained spectrum.  It is whether an independent axiom
selects one projection, or whether an analytic universality theorem can be
proved despite the negative finite flow.

## Status ledger

- **DERIVED:** all three f-vectors and all exact recovery sums.
- **DERIVED:** exact repeated maximum local differences `1/15` and `5/32`.
- **DERIVED NEGATIVE:** no coefficientwise supremum convergence from `k=2`
  to `k=4`.
- **DERIVED NUMERICAL:** all strong-limit Betti nullities and Ritz residuals.
- **DERIVED ON CONTROLS:** degrees 0 and 1 have identical constrained pencils.
- **PATTERN:** degree-2 first-quartet isospectrality despite different
  matrices.
- **PATTERN NEGATIVE:** top-degree relative separation grows by factor 13.72.
- **OPEN:** common or distinct analytic continuum limits.
- **OPEN:** a physical principle selecting `P^C`, `P^D`, or neither.
- **NOT CLAIMED:** round-`S^3` geometry, Lorentzian time, causal speed,
  inertia, mass or Planck units.

## Reproduction

```bash
/home/razvan/science/.venv/bin/python -u \
  reproducible/verify_whitney_broken_feec_universality.py
```

Expected result: `8/8`.

