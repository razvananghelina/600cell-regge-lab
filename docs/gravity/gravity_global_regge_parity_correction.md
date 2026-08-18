# Disclosed correction: ordered phase parity changes the regular Hessian

Date: 2026-08-12

Original orbit-reduction protocol commit: `06a1c6a`

Status: **correction frozen after the first complete run returned 38/39**.

## Observed result

Every preregistered construction and consistency control passed.  The only
failure was a verifier assertion, stronger than the protocol required, that
the two ordered phase-parity classes should have equal regular Hessian
singular spectra.

The observed relative maximum spectral separation was

```text
1.6233464461549153e-2.
```

The numerical controls at the same step were:

```text
even Hessian symmetry residual  3.11e-10
odd  Hessian symmetry residual  3.06e-10
maximum imaginary residual      3.99e-9
rank at 1e-7,1e-9,1e-11         35 for both
even s_min/s_max                 1.8081e-2
odd  s_min/s_max                 1.5653e-2.
```

The actions and gradient norms at regular data agree between parities to
displayed precision, so the distinction first appears at second order.

The original written protocol asked to determine whether parity changes the
linearized algebra; it did not preregister equality as the acceptance answer.
The implementation nevertheless encoded equality as a PASS condition.  That
was a framing error.  It will not be silently changed: this correction records
the seen rank and separation before the condition is replaced by a decisive
separation test.

## Frozen independent controls

The corrected verifier must retain all original 38 passing checks and add:

1. report signed Hessian eigenvalues as well as singular values for both
   parities;
2. report basis-invariant trace, Frobenius norm and signed log-determinant;
3. compute the Hessian trace independently from centered second differences
   of the complete reduced action in all 35 coordinate directions at step
   `5e-4`;
4. require that direct trace to agree with the gradient-derived Hessian trace
   to relative tolerance `2e-5` for both parities;
5. call parity equality only if the relative singular-spectrum difference is
   below `2e-5`, and call parity separation only if it exceeds `1e-3` while
   the controls above pass;
6. record the original `38/39` run and this correction protocol hash in the
   result JSON and note.

The interval between `2e-5` and `1e-3` is deliberately inconclusive.  The
observed `1.62e-2` lies well outside it, but the value must be recomputed.

## Interpretation boundary

If confirmed:

- **DERIVED COMPUTATIONAL:** the two staircase phase parities have different
  linearized Regge actions at the common regular metric;
- **DERIVED NEGATIVE:** phase order is not erased by the present discrete
  Regge action at quadratic order;
- **OPEN:** whether a continuum/refinement limit restores refoliation
  invariance;
- **OPEN:** whether either parity is dynamically selected;
- **NOT CLAIMED:** that the parity is fermion chirality, CP, or the arrow of
  time.

This correction still performs no stationary-root search.
