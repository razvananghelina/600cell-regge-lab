# Preregistration: does the existing `Tr(Box^3)` constraint select six fibrations?

Date: 2026-08-10

This protocol is committed before checking stationarity or computing the
critical/level set.  It follows the open boundary in
`hopf_box_projector_lift_verdict.md`.

One failed exploratory constrained optimizer produced overflow and no usable
result.  No numerical extremum or critical value from that run is accepted as
evidence or used below.

## Fixed data

Use exactly the five-dimensional integer operator space

```text
W = span_R{Box_F : F is one of the six certified Hopf fibrations},
Box_F = 6 A_f,F - A.
```

The already certified identities are

```text
sum_F Box_F = 0,
Tr(Box_F^2) = 7200,
Tr(Box_F^3) = 14400.
```

No additional invariant, weight, or normalization may be introduced.

Define on `W`

```text
q(X) = Tr(X^2),
f(X) = Tr(X^3).
```

## Ordered tests

1. Test exactly whether each `Box_F` is stationary for `f` on the sphere
   `q=7200`.  This means checking the full differential equation

   ```text
   df_X = lambda dq_X  on all of W,
   ```

   not merely along the other five simplex vertices.
2. If a `Box_F` is not stationary, exhibit an exact tangent direction with
   `dq=0` and `df!=0`.  Then use the regular level-set/implicit-function
   criterion to conclude that `q=7200, f=14400` contains a local real
   continuum through that vertex.  This immediately refutes selection.
3. If all six vertices are stationary, compute the complete exact critical
   scheme of `f|_{q=7200}`, its real critical values and orbit sizes.  Do not
   infer uniqueness from numerical optimization.
4. Independently determine whether the simultaneous level set
   `q=7200, f=14400` contains anything besides the six vertices.  An explicit
   exact extra point is sufficient to refute uniqueness; a positive result
   requires an exhaustive algebraic proof.

## Decision boundary

- **Advance:** the already derived constraints select exactly the six
  `Box_F`, with an exhaustive exact proof.
- **Kill:** any additional point, curve, or higher-dimensional component on
  the same level set.  In particular, nonstationarity of a desired vertex is
  already a local-continuum kill.

Even an advance would remain an operator-selection theorem, not yet a
licensed spectral-action or physical vacuum theorem.
