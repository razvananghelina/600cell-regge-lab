# Protocol: coefficient-ball root count for the negative-shape recurrence

Date: 2026-08-18

Prior-art gate commit: `3467728`

This protocol is **TARGET-DISCLOSED**.  The preceding committed calculation
already reported the midpoint pattern `15` roots inside and `15` outside the
unit circle.  The present test is not evidence for that midpoint pattern.  Its
only new question is whether a complete coefficient-ball certificate transfers
a root count to every admissible recurrence.

No root count produced after this commit may alter the carrier, error balls,
contour, safety thresholds, subdivision schedule or outcome ordering below.

## 1. Frozen object and complete hypothesis list

For each of

```text
parity  = even, odd
sector  = 4, 5
variant = operational_primary, operational_shadow,
          validation_primary, validation_shadow,
```

reconstruct exactly the `15`-dimensional negative-stiffness carrier selected by
the blind stiffness census and inherited by the autonomous-dynamics verifier.
Use the same centered archive, conformal/shape restrictions, eigengap bounds,
and invariance bounds.  No basis or tolerance is recomputed from the root
spectrum.

On each carrier audit

```text
Q(z) = A2 z^2 + A1 z + A0,
A2   = I + Gamma,
A1   = -2 I + Omega,
A0   = I - Gamma.
```

The hypotheses are:

1. all committed source hashes and upstream outcomes match;
2. the conformal, shape and negative carriers remain resolved;
3. both `Gamma` and `Omega` leave the negative carrier
   invariant-consistently;
4. the inherited operator-norm balls are
   `||DeltaGamma|| <= epsilon_Gamma` and
   `||DeltaOmega|| <= epsilon_Omega`, with the residual-to-carrier terms added
   exactly as in the upstream verifier;
5. `I+Gamma` is regular-resolved throughout its ball, so the polynomial has
   exactly `30` finite roots counting algebraic multiplicity.

This is a local frozen-coefficient diagnostic.  Time stationarity, reciprocal
pairing, later-slab evolution, nonlinear stability and refinement are not
hypotheses.

## 2. Continuous Rouche bound

On `z=exp(i theta)`, preserve the correlated perturbation

```text
DeltaQ(z) = DeltaGamma (z^2 - 1) + DeltaOmega z.
```

Thus

```text
delta(theta) = epsilon_Gamma |z^2-1| + epsilon_Omega.
```

Let

```text
evaluation_floor = 1000 eps_machine m
                   max(1, 2||A2|| + ||A1|| + ||A0||),
L_Q = 2||A2|| + ||A1||,
L_delta = 2 epsilon_Gamma,
```

where `m=15` and all norms are spectral norms.  For an angular interval with
centre `theta_c` and half-width `h`, Weyl's singular-value inequality and the
derivative bounds give

```text
min_interval [sigma_min(Q)-s(delta+evaluation_floor)]
 >= sigma_min(Q(theta_c))
    - s(delta(theta_c)+evaluation_floor)
    - (L_Q+s L_delta) h.
```

The deterministic cover begins with `256` equal closed angular intervals,
bisects only unresolved intervals, uses depth-first left-before-right order,
and stops at maximum subdivision depth `32` or `2,000,000` evaluated intervals
per cell and safety.  If the pointwise midpoint margin itself is nonpositive,
the cover fails immediately because subdivision cannot repair that witness.

Run the cover twice:

- `s=1`, the literal matrix-valued Rouche condition;
- `s=100`, the repository's resolved-vs-error safety standard.

No angular sample may be deleted.  A failed cover is `OPEN`; it is not evidence
that a root lies on the circle.

## 3. Independent midpoint counts

The midpoint count is not transferred unless both methods below agree.

### 3.1 Generalized/companion count

Use the inherited regular `A2` to construct

```text
C = [[0, I], [-A2^(-1) A0, -A2^(-1) A1]].
```

Count its `30` eigenvalues strictly inside and outside the unit circle.  Record
the maximum normalized polynomial residual of every computed eigenpair.  It
must be below `1e-10`.

### 3.2 Argument-principle count

Evaluate the complex sign of `det Q(exp(i theta))` on uniform closed meshes of
`8192` and `16384` intervals.  Unwrap its phase.  Both meshes must yield the
same integer winding, the winding must equal the companion inside count, and
the largest principal phase increment on the finer mesh must be below `pi/4`.
The continuous Rouche cover separately proves that the determinant never
vanishes on the contour; this winding computation is only the independent
midpoint count.

## 4. Controls

Before any inherited cell is accepted:

1. `Q_free(z)=(z-1)^2 I_2` must fail the literal unit-circle cover;
2. `Q_unit(z)=(z^2+1) I_2` must fail the literal cover;
3. `Q_hyper(z)=(z^2-2.5z+1) I_2`, with coefficient errors `1e-12`, must pass
   the `s=100` cover and be counted as two inside and two outside roots;
4. reversing the angular enumeration must leave every synthetic winding count
   unchanged after correcting contour orientation.

## 5. Preregistered ledger and outcomes

For all `16` inherited cells write:

- coefficient norms and errors;
- leading-coefficient regularity margin;
- literal and `100x` cover status, interval count, maximum depth, weakest lower
  margin, and any first pointwise failure;
- both midpoint counts, phase-resolution diagnostics and polynomial residual;
- the transferred counts inside/on/outside if and only if the cover passes.

Outcome order:

1. `NEGATIVE_SHAPE_ROOT_COUNT_CONTROL_FAILED` if a hash, carrier, regularity or
   synthetic control fails;
2. `NEGATIVE_SHAPE_ROOT_COUNT_MIDPOINT_INCONSISTENT` if the two midpoint counts
   disagree, a residual gate fails, or a phase gate fails;
3. `NEGATIVE_SHAPE_ROOT_COUNT_COEFFICIENT_OPEN` if any literal (`s=1`) cover
   fails;
4. `NEGATIVE_SHAPE_ROOT_COUNT_SAFETY_OPEN` if all literal covers pass but any
   `s=100` cover fails;
5. `NEGATIVE_SHAPE_ROOT_COUNT_SCHEDULE_DEPENDENT` if all covers pass but the
   transferred counts differ among derivative variants or parities for a
   sector;
6. `NEGATIVE_SHAPE_LOCAL_HYPERBOLIC_RESOLVED` if all `s=100` covers pass and
   every cell has at least one root strictly inside and at least one strictly
   outside the unit circle;
7. `NEGATIVE_SHAPE_LOCAL_NONEXPANDING_RESOLVED` only if all `s=100` covers pass
   and every cell has zero roots outside;
8. `NEGATIVE_SHAPE_ROOT_COUNT_OTHER_RESOLVED` for any other complete,
   schedule-robust transferred count.

The word `LOCAL` is mandatory in any positive hyperbolicity interpretation.
Nothing here establishes a cosmological arrow, a physical growth rate, a
long-time instability, a continuum gravitational wave or a value of `c`.

