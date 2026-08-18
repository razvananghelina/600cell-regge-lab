# Result: the finite-family generalized bundle rotates between ticks

Date: 2026-08-18

## Headline

The high-precision a posteriori audit resolves a nonzero rotation of the
rank-15 generalized-mode fiber between the old and shifted centered ticks in
both disclosed symmetry sectors and both schedule parities.

All `64 = 2 parities * 2 sectors * 4 old schedules * 4 shifted schedules`
old/new comparisons are

```text
ROTATION_RESOLVED.
```

The old/new Euclidean projector distance is

```text
3.960493055435e-7,
```

while the complete projector comparison error is approximately

```text
7.20e-53 ... 8.72e-53.
```

The complete within-time diameter of the four derivative schedules is only

```text
4.80e-25 ... 4.90e-25.
```

Thus the old/new separation is approximately

```text
8.09e17 ... 8.26e17
```

times the complete finite-family diameter, far beyond the preregistered factor
`100`.  The frozen outcome is

```text
RESIDUAL_FINITE_FAMILY_ROTATION_RESOLVED.
```

This refutes the literal reading that the two finite-family fibers are the
same subspace.  It does not yet derive a physical connection or transport
between the distinct fibers.

## Provenance ledger

| stage | commit |
|---|---|
| primary-literature and framing gate | `7e815de` |
| original preregistration | `1b48e52` |
| pre-execution finite-family robustness correction | `ad2eadc` |
| registered verifier before first execution | `c69a783` |
| serializer-only repair after failed artifact write | `8f9e755` |
| deterministic result artifact | `61eb381` |

The result artifact has SHA-256

```text
3244185127aecf7c9a44261cced0be521c9dc42bf8e44f909d8a0ce10a96eadf.
```

The accepted verifier source has SHA-256

```text
ccf2ebe03c6e39c3d6e6b538d1c02d278804553987d65db0eeb67fce7936ca5a.
```

The first scientific execution completed the numerical calculation and
disclosed the preregistered outcome, but failed at the final JSON write because
a `numpy.bool_` was not serializable.  The repair cast only affected artifact
fields to built-in `bool`; it changed no scientific value or branch.  The
complete calculation was then executed twice, passed `10/10` both times and
produced the byte-identical artifact above.  The failed execution is not
counted as a reproducible result.

Only this targeted verifier was run.  The full suite was not run.

## What changed relative to the broad audit

The previous audit collapsed all four derivative schedules into one
unstructured whole-pencil error ball.  That allowed its entire uncertainty to
occupy the off-diagonal spectral block, producing projector errors of
approximately `2.43e-4 ... 6.46e-4`.  Those errors were much larger than the
observed old/new displacement, so “common fiber” meant only zero-consistency.

The present audit does not choose a smaller tolerance.  It:

1. retains all four disclosed finite-difference schedules separately;
2. reconstructs their Flint source balls before binary64 conversion;
3. builds the rank-5 conformal image and rank-25 action-selected shape
   nullspace at 100 decimal digits;
4. Cholesky-whitens every Hermitian-definite pencil;
5. certifies the computed lower spectral subspace from its actual
   off-diagonal residual and gap;
6. compares every old schedule with every shifted schedule;
7. requires the minimum old/new separation to exceed the full within-family
   diameter by a factor greater than `100`.

No eigenvector matching, Procrustes alignment, polar transport or fitted basis
was used.  Projectors make the result independent of rotations inside the
rank-15 cluster.

## Numerical certificate

Every one of the 32 schedule-local cells retains:

```text
shape dimension              25,
generalized lower cluster    15,
generalized upper cluster    10,
positive restricted kinetic form,
generalized gap              3.997703589e-5 ... 3.997706876e-5.
```

The relevant error ranges are:

| quantity | range |
|---|---:|
| Flint `M` source radius | `5.29e-58 ... 9.48e-57` |
| Flint `V` source radius | `1.17e-57 ... 1.93e-56` |
| whitened off-diagonal residual | `3.80e-105 ... 6.08e-105` |
| schedule-local projector error | `7.61e-54 ... 7.85e-53` |
| within-time schedule distance | `2.62e-35 ... 4.90e-25` |
| old/new schedule distance | `3.960493055435e-7` |
| old/new distance/error units | `4.54e45 ... 5.50e45` |

The high-precision projectors also overlap all previously committed binary64
projector controls; their binary distances are `1.56e-14 ... 2.36e-13`.

All 48 within-time schedule pairs are themselves numerically distinguishable
at this precision, as finite-step approximations should be.  That fact is not
called physical scheme dependence: their complete diameter is 18 orders of
magnitude smaller than the old/new separation.  The preregistered evidence is
the scale separation, not exact equality of finite-step formulas.

## Post-result literature check

The post-result terminology leads directly to Grubišić, Truhar and Veselić,
*The rotation of eigenspaces of perturbed matrix pairs*, which develops
`sin(Theta)` bounds for positive-definite matrix pairs and emphasizes that
sharpness depends on the matrix-dependent scalar product:
<https://doi.org/10.1016/j.laa.2012.01.026> and
<https://arxiv.org/abs/1011.4424>.

This does not supply the present 600-cell result.  It confirms that the
generalized-pair rotation is a standard mathematical object and sharpens one
caveat: the intermediate eigenspace estimate naturally lives in the kinetic
metric, whereas the reported final distance is the ordinary Euclidean
projector distance after the explicit Cholesky and shape lift.  The metric
conversion is retained in the error bound rather than silently identifying
the two angles.

The external novelty of applying this construction to the present Regge--dust
trajectory remains **OPEN** pending a dedicated review.

## Status ledger

- **DERIVED COMPUTATIONAL, conditional on the frozen four-schedule derivative
  family:** the old and shifted rank-15 generalized fibers are distinct in all
  64 cross-schedule comparisons.
- **DERIVED COMPUTATIONAL:** the old/new displacement exceeds the full
  structured finite-family diameter by approximately `8e17` in all four
  parity/sector cells.
- **DERIVED COMPUTATIONAL:** the earlier common-fiber label was caused by an
  intentionally conservative unstructured family error, not by an equal
  high-precision midpoint projector.
- **STRUCTURAL:** the rank-30 object should now be treated as a time-dependent
  bundle consisting of two rank-15 symmetry sectors, not as one fixed
  configuration subspace.
- **OPEN:** an analytic/automatic enclosure of the exact Hessian outside the
  four finite-step schedules.
- **OPEN:** an action-selected connection or transport between the old and
  shifted fibers.
- **OPEN:** whether the full Jacobi recurrence maps the old fiber into the
  shifted fiber.
- **OPEN:** reduced propagator, characteristic roots, long-time stability,
  continuum/refinement behavior, gravitons, waves, particle inertia, mass and
  limiting speed.

Here “inertia” still means the signature of the quadratic Regge action.  This
result is not a derivation of ordinary particle inertia.

## Next load-bearing gate

Because the fibers are distinct, the literal identity is not a valid
trivialization for a reduced two-step product.  A polar or Procrustes alignment
would be mathematically natural but physically fitted: the Regge action has
not selected it.

The next test must therefore use the already derived full tangent/Jacobi map
itself.  Lift the old rank-15 fiber to the appropriate 30-position/30-momentum
phase-space data, apply the action-generated non-autonomous tangent step, and
measure its component outside the shifted rank-15 fiber.  Only if that leakage
is zero-consistent may the restricted map be called an action-selected bundle
transport and used to form a reduced propagator.
