# Result: literal root count transfers, but the project safety margin remains open

Date: 2026-08-18

## Headline

**DERIVED COMPUTATIONAL, TARGET-DISCLOSED.**  For every one of the `16`
inherited parity/sector/derivative cells, the continuous matrix-valued Rouche
test succeeds with the literal coefficient balls.  Both independent midpoint
counts give

```text
inside |z|<1    15
on     |z|=1     0
outside|z|>1    15
```

and Rouche transfers that count to every polynomial inside each declared ball.

The repository's deliberately stronger `100 x error` rule fails in all `16`
cells.  The preregistered outcome is therefore

```text
NEGATIVE_SHAPE_ROOT_COUNT_SAFETY_OPEN,
```

not `NEGATIVE_SHAPE_LOCAL_HYPERBOLIC_RESOLVED`.

This is not a contradiction.  The literal theorem requires a factor strictly
larger than `1`; the smallest sampled ratio is approximately `8.16755`.  That
is enough for the stated coefficient balls, but not enough for the project's
extra factor `100` against unmodelled numerical error.

## Provenance

| stage | commit |
|---|---|
| primary-literature and structural gate | `3467728` |
| target-disclosed protocol | `925811e` |
| verifier registered before first execution | `a9d3796` |
| pre-spectrum sparse-import fix | `d07443c` |
| pre-spectrum multiprecision-unit fix | `a4873db` |
| post-result margin-ledger clarification | `b231cf8` |

The first two executions stopped before constructing any target polynomial:
one extracted audited helper required the omitted global `scipy.sparse`, and
the high-precision sector helper required the omitted multiprecision unit `I`.
Both failures and fixes precede the first root count.

Verifier:

```text
reproducible/verify_gravity_600cell_dust_negative_shape_root_count.py
```

Artifact:

```text
reproducible/gravity_600cell_dust_negative_shape_root_count.json
SHA-256 7f71d680f4ba34da2f3a8af7c8e4b92668e1cf16e682acbf83f9d6bbec41b280
```

Two final targeted runs were byte-identical and each reported `13/13`.  The
full suite was deliberately not run, following the user's instruction.

## What was actually certified

On each negative carrier the frozen recurrence has polynomial

```text
Q(z) = (I+Gamma) z^2 + (-2I+Omega) z + (I-Gamma).
```

For `z=exp(i theta)`, the inherited correlated perturbation is

```text
DeltaQ(z) = DeltaGamma (z^2-1) + DeltaOmega z.
```

The adaptive deterministic cover proves throughout the unit circle that

```text
sigma_min(Q(z))
  > epsilon_Gamma |z^2-1| + epsilon_Omega + evaluation_floor.
```

The Gohberg--Sigal/Rouche theorem then preserves the enclosed algebraic root
count for every admitted perturbation.  Regularity of `I+Gamma` proves that
there are `30` finite roots in total; hence `15` enclosed roots imply `15`
outside roots and no unit-circle roots.

All literal covers required exactly `18,368` evaluated angular intervals and
reached maximum depth `16`.  The weakest certified leaf lower bound lies
between `2.8112467e-10` and `2.8112509e-10`.  This tiny leaf number is an
artifact of stopping subdivision immediately after a positive bound; the more
informative smallest sampled signal/error ratio is

```text
8.1675481076 ... 8.1675481114.
```

The inherited errors are nearly identical across all cells:

```text
epsilon_Gamma approximately 4.88735e-8
epsilon_Omega approximately 1.41262e-7.
```

The `100x` cover encounters its first pointwise failure at approximately
`theta=0.00306796`, with ratio `74.632`.  That is only the first failure of its
left-to-right traversal, not the global minimum; the completed literal cover
later samples the smaller ratio `8.17`.

## Independent midpoint checks

Every cell gives the same count by both methods:

1. the regular quadratic companion has `15/0/15` roots
   inside/on/outside;
2. the winding of `det Q` on both the `8192`- and `16384`-interval meshes is
   `15`.

The largest normalized polynomial residual is below `9.05e-16`.  The largest
principal phase increment is about `0.0191`, far below the preregistered
`pi/4` gate.  The nearest midpoint root is approximately `0.00107064` in
modulus away from the unit circle.

The free double-unit-root control and the `z^2+1` unit-root control both reject
the contour certificate.  The scalar hyperbolic control with roots `0.5` and
`2` passes at `100x` and its forward/reverse windings agree after orientation
correction.

## Hostile interpretation audit

1. **LOCAL ONLY.**  This is the root count of a frozen two-step recurrence.
   It is not a long-time Lyapunov exponent and does not show that independently
   solved future slabs repeat the same coefficients.
2. **AUTONOMY IS CALIBRATED, NOT EXACT.**  The two negative spaces are
   invariant-consistent within the inherited error model.  No exact symbolic
   theorem yet says that the nonlinear theory preserves them.  Coupling to the
   complementary shape modes remains physically **OPEN**.
3. **THE `100x` FAILURE MATTERS.**  The result is mathematically stable inside
   the declared balls, but only by a worst sampled factor `8.17`.  It would be
   misleading to call this a resolved physical instability under the project's
   established numerical semantics.
4. **NO RECIPROCAL-PAIR THEOREM WAS USED.**  The inherited time-asymmetry
   coefficient `N` is Hermitian-consistent, not anti-Hermitian; the full-source
   relative palindromic defect is approximately `0.334`.  The observed equal
   inside/outside count is an output, not a structural assumption.
5. **NO POLARIZATION HAS BEEN IDENTIFIED.**  The `30` real position modes have
   not been classified as scalar, vector or tensor perturbations, and no
   continuum gravitational-wave dispersion or speed follows.

## Status ledger

- **DERIVED COMPUTATIONAL:** the literal coefficient-ball Rouche condition
  holds in all `16` cells.
- **DERIVED COMPUTATIONAL:** every declared ball has root count `15/0/15`.
- **DERIVED UPSTREAM:** the carrier contains `30` real negative-stiffness
  positions and is invariant-consistent under the frozen dynamics.
- **OPEN BY PREREGISTERED SAFETY RULE:** the minimum observed margin factor is
  `8.17`, below `100`.
- **OPEN:** exact nonlinear invariance, later-slab persistence, Lyapunov growth,
  refinement, continuum polarization and physical speed.
- **REFUTED AS A SHORTCUT:** ordinary companion Bauer--Fike bounds and the
  companion singular norm do not decide the root count; direct polynomial
  structure is materially stronger.

## Next load-bearing gate

There are two defensible continuations, neither of which may be chosen by
looking for a desired count:

1. tighten the source derivative and carrier balls by an independently
   preregistered higher-precision reconstruction; or
2. use the already selected kinetic form as a canonical metric and derive a
   mass-weighted Rouche bound.  This is admissible only if the metric is fixed
   before the spectrum and if its condition-number transport is included.

An arbitrary similarity transformation optimized to enlarge `8.17` would be
fitting and is forbidden.  Even if the `100x` gate is later crossed, the next
physical question is non-autonomous evolution across independently solved
slabs, not another frozen diagonalization.

